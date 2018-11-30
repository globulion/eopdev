//#include "psi4/libtrans/integraltransform.h"
//#include "psi4/libdpd/dpd.h"

#include "solver.h"

using namespace std;
using namespace psi;
using namespace oepdev;

using SharedMTPConv = std::shared_ptr<oepdev::MultipoleConvergence>;

// ===> Repulsion Energy <=== //
RepulsionEnergySolver::RepulsionEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Benchmarks
  methods_benchmark_.push_back("HAYES_STONE"     );
  methods_benchmark_.push_back("DENSITY_BASED"   );
  methods_benchmark_.push_back("MURRELL_ETAL"    );
  methods_benchmark_.push_back("OTTO_LADIK"      );
  methods_benchmark_.push_back("EFP2"            );
  // OEP-based
  methods_oepBased_ .push_back("MURRELL_ETAL_GDF_ESP");
  methods_oepBased_ .push_back("MURRELL_ETAL_GDF_CAMM");
  methods_oepBased_ .push_back("MURRELL_ETAL_ESP");
}
RepulsionEnergySolver::~RepulsionEnergySolver() {}
double RepulsionEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "MURRELL_ETAL_GDF_ESP")  e = compute_oep_based_murrell_etal_gdf_esp();
  else if (method == "MURRELL_ETAL_GDF_CAMM") e = compute_oep_based_murrell_etal_gdf_camm();
  else if (method == "MURRELL_ETAL_ESP"    )  e = compute_oep_based_murrell_etal_esp();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for repulsion energy calculations!\n");
  }
  return e;
}
double RepulsionEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" ||
           method == "HAYES_STONE"  ) e = compute_benchmark_hayes_stone();
  else if (method == "DENSITY_BASED") e = compute_benchmark_density_based();
  else if (method == "MURRELL_ETAL" ) e = compute_benchmark_murrell_etal();
  else if (method == "OTTO_LADIK"   ) e = compute_benchmark_otto_ladik();
  else if (method == "EFP2"         ) e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
double RepulsionEnergySolver::compute_benchmark_density_based() {
  double e             = 0.0;
  double e_Pauli_nuc   = 0.0;
  double e_Pauli_Pauli = 0.0;
  double e_Pauli_el    = 0.0;
  double e_Pauli_kin   = 0.0;
  double e_exch        = 0.0;

  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Density-Based (2012))");
  psi::timer_on("Solver E(Paul) Density-Based    ");

  // ---> Allocate <--- //
  int nbf   = wfn_union_->basisset()->nbf();
  psi::IntegralFactory fact(wfn_union_->basisset());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> kinInt(fact.ao_kinetic());
  std::shared_ptr<psi::OneBodyAOInt> potInt(fact.ao_potential());

  std::shared_ptr<psi::Matrix> Sao = std::make_shared<psi::Matrix>("Sao", nbf, nbf);
  std::shared_ptr<psi::Matrix> Tao = std::make_shared<psi::Matrix>("Tao", nbf, nbf);
  std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao", nbf, nbf);

  // ---> Compute One-Electron Integrals <--- //
  ovlInt->compute(Sao);
  kinInt->compute(Tao);
  potInt->compute(Vao);

  // ---> Compute Density Matrix after Orthogonalization <--- //
  std::shared_ptr<psi::Matrix> Ca_occ = wfn_union_->Ca_subset("AO", "OCC");
  std::shared_ptr<psi::Matrix> Smo = psi::Matrix::triplet(Ca_occ, Sao, Ca_occ, true, false, false);
  Smo->invert();
  std::shared_ptr<psi::Matrix> Da_ao_oo = psi::Matrix::triplet(Ca_occ, Smo, Ca_occ, false, false, true);

  // ---> Compute Difference Pauli Density Matrix <--- //
  std::shared_ptr<psi::Matrix> Imo = Smo->clone();
  Imo->identity();
  Smo->subtract(Imo);
  std::shared_ptr<psi::Matrix> DeltaDa_ao = psi::Matrix::triplet(Ca_occ, Smo, Ca_occ, false, false, true);

  // ---> One-Electron Contribution <--- //
  e_Pauli_nuc = 2.0 * DeltaDa_ao->vector_dot(Vao);
  e_Pauli_kin = 2.0 * DeltaDa_ao->vector_dot(Tao);

  // ---> Two-Electron Contribution <--- //
  std::shared_ptr<psi::IntegralFactory> ints = std::make_shared<psi::IntegralFactory>
                                               (wfn_union_->basisset(), wfn_union_->basisset(),
                                                wfn_union_->basisset(), wfn_union_->basisset());

  std::shared_ptr<psi::TwoBodyAOInt> tei(ints->eri());
  const double * buffer = tei->buffer();
  std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(ints, "ALL");
  int i, j, k, l; int nbf_1 = wfn_union_->l_nbf(0); int nbf_m = nbf_1-1;
  double integral;
  double** dD  = DeltaDa_ao->pointer();
  double** da  = Da_ao_oo->pointer();
  double** Da  = wfn_union_->Da()->pointer();

  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       shellIter->compute_shell(tei);
       std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            i = intsIter->i();
            j = intsIter->j();
            k = intsIter->k();
            l = intsIter->l();

            integral = buffer[intsIter->index()];
            e_Pauli_el    += dD[i][j] * Da[k][l] * integral;
            e_Pauli_Pauli += dD[i][j] * dD[k][l] * integral;
            e_exch        -= da[i][l] * da[j][k] * integral;
            if (((i < nbf_1) && (j < nbf_1) && (k < nbf_1) && (l < nbf_1)) ||
                ((i > nbf_m) && (j > nbf_m) && (k > nbf_m) && (l > nbf_m)))
                  e_exch  += Da[i][l] * Da[j][k] * integral;
       }
  }
  e_Pauli_el    *= 4.0;
  e_Pauli_Pauli *= 2.0;

  // ---> Finish <--- //
  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Density-Based (2012))");
  psi::timer_off("Solver E(Paul) Density-Based    ");
  e = e_Pauli_nuc + e_Pauli_kin + e_Pauli_el + e_Pauli_Pauli;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (Density-Based)              <==\n\n");
     psi::outfile->Printf("     PAU NUC   = %13.6f\n", e_Pauli_nuc                      );
     psi::outfile->Printf("     PAU KIN   = %13.6f\n", e_Pauli_kin                      );
     psi::outfile->Printf("     PAU EL    = %13.6f\n", e_Pauli_el                       );
     psi::outfile->Printf("     PAU PAU   = %13.6f\n", e_Pauli_Pauli                    );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP_1   = %13.6f\n", e_Pauli_nuc+e_Pauli_kin          );
     psi::outfile->Printf("     E REP_2   = %13.6f\n", e_Pauli_Pauli+e_Pauli_el         );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch                           );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch                         );
     psi::outfile->Printf("\n");
  }

  // ---> Return Total Pauli Repulsion Energy <--- //
  return e;
}
double RepulsionEnergySolver::compute_benchmark_hayes_stone() {

  // ===> Start computations <=== //
  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Hayes-Stone (1984))");
  psi::timer_on("Solver E(Paul) Hayes-Stone      ");

  double e_1, e_2, e_ex;

  // ===> One electron part <=== //

  int nbf   = wfn_union_->basisset()->nbf();

  // ---> Allocate <--- //
  psi::IntegralFactory fact(wfn_union_->basisset());

  std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> kinInt(fact.ao_kinetic());
  std::shared_ptr<psi::OneBodyAOInt> potInt(fact.ao_potential());

  std::shared_ptr<psi::Matrix> Sao       = std::make_shared<psi::Matrix>("Sao"  , nbf, nbf);
  std::shared_ptr<psi::Matrix> Tao       = std::make_shared<psi::Matrix>("Tao"  , nbf, nbf);
  std::shared_ptr<psi::Matrix> Vao       = std::make_shared<psi::Matrix>("Vao"  , nbf, nbf);

  // ---> Compute one-electron integrals <--- //
  ovlInt->compute(Sao);
  kinInt->compute(Tao);
  potInt->compute(Vao);

  // ---> Accumulate one electron part <--- //
  Tao->add(Vao);

  // ---> Transform one electron and overlap matrices to MO basis <--- //
  std::shared_ptr<psi::Matrix> Ca_occ = wfn_union_->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Smo = psi::Matrix::triplet(Ca_occ, Sao, Ca_occ, true, false, false);
  std::shared_ptr<psi::Matrix> Tmo = psi::Matrix::triplet(Ca_occ, Tao, Ca_occ, true, false, false);

  // ---> Invert the overlap matrix in MO basis <--- //
  Smo->invert();

  // ---> Finalize with one-electron term <--- //  
  std::shared_ptr<psi::Matrix> Imo = std::make_shared<psi::Matrix>(Smo);
  Imo->identity();
  Smo->subtract(Imo);
  e_1 = 2.0 * Tmo->vector_dot(Smo);

  // ===> Two electron part <=== //
  Smo->add(Imo);
  // ---> Loop over ERI's in MO space <--- //
  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  dpdbuf4 buf;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[O,O]"  ), integrals->DPD_ID("[O,O]"  ),
                          integrals->DPD_ID("[O>=O]+"), integrals->DPD_ID("[O>=O]+"), 0, "MO Ints (OO|OO)");

  double integral, v;
  double** T = Smo->pointer();
  e_2 = 0.0;
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
       for (int kl = 0; kl < buf.params->rowtot[h]; ++kl) {
            int k = buf.params->roworb[h][kl][0];
            int l = buf.params->roworb[h][kl][1];
            for (int mn = 0; mn < buf.params->coltot[h]; ++mn) {
                 int m = buf.params->colorb[h][mn][0];
                 int n = buf.params->colorb[h][mn][1];
                 //psi::outfile->Printf(" Jint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);

                 integral = buf.matrix[h][kl][mn];

                                           v = 2.0*T[k][l]*T[m][n] - T[k][n]*T[l][m];
               //if  (k != l)              v+= 2.0*T[k][l]*T[m][n] - T[l][n]*T[k][m]; 
               //if  (m != n)              v+= 2.0*T[k][l]*T[m][n] - T[l][n]*T[k][m];
               //if ((k != l) && (m != n)) v+= 2.0*T[k][l]*T[m][n] - T[k][n]*T[l][m];
                 if ((k == l) && (m == n)) v-= 2.0;
                 if ((k == n) && (l == m)) v+= 1.0;

                 e_2 += integral * v;
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);

  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Hayes-Stone (1984))");
  psi::timer_off("Solver E(Paul) Hayes-Stone      ");

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  // ===> Compute Exchange Energy <=== //
  e_ex = compute_pure_exchange_energy();

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (Hayes-Stone)                <==\n\n");
     psi::outfile->Printf("     E REP 1   = %13.6f\n", e_1                              );
     psi::outfile->Printf("     E REP 2   = %13.6f\n", e_2                              );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e_1+e_2                          );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_ex                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e_1+e_2+e_ex                     );
     psi::outfile->Printf("\n");
  }

  // ---> Return Total Repulsion Energy <--- //
  return e_1+e_2+e_ex;
}
double RepulsionEnergySolver::compute_pure_exchange_energy() {
  //psi::timer_on ("SOLVER: HF Exchange Energy Calculations");
  psi::timer_on("Solver E(Exch) Hayes-Stone      ");

  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  dpdbuf4 buf_1212;

  global_dpd_->buf4_init(&buf_1212, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[1,2]"  ),
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[1,2]"  ), 0, "MO Ints (12|12)");

  double e_ex = 0.0;
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1212, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1212, h);
       for (int pq = 0; pq < buf_1212.params->rowtot[h]; ++pq) {
            int p = buf_1212.params->roworb[h][pq][0];
            int q = buf_1212.params->roworb[h][pq][1];
            for (int rs = 0; rs < buf_1212.params->coltot[h]; ++rs) {
                 int r = buf_1212.params->colorb[h][rs][0];
                 int s = buf_1212.params->colorb[h][rs][1];
                 if ((p==r) && (q==s)) e_ex += buf_1212.matrix[h][pq][rs];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1212.matrix[h][pq][rs]);
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1212, h);
  }
  global_dpd_->buf4_close(&buf_1212);
  e_ex *= -2.0;

  //psi::timer_off("SOLVER: HF Exchange Energy Calculations");
  psi::timer_off("Solver E(Exch) Hayes-Stone      ");

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  return e_ex;
}
double RepulsionEnergySolver::compute_benchmark_murrell_etal() {
  double e_s1 = 0.0, e_s2 = 0.0, e = 0;

  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell et al.)");
  psi::timer_on("Solver E(Paul) Murrell-etal     ");
  int nbf_1 = wfn_union_->l_nbf(0);
  int nbf_2 = wfn_union_->l_nbf(1);

  // ===> One electron part <=== //
  std::shared_ptr<psi::Matrix> VaoA21    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_1);
  std::shared_ptr<psi::Matrix> VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> VaoA22    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_2);
  std::shared_ptr<psi::Matrix> VaoB11    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)"  , nbf_1, nbf_2);
  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_11(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_22(wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));
  std::shared_ptr<psi::OneBodyAOInt> oneInt, ovlInt(fact_12.ao_overlap());
  std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_12.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_21.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_11= std::make_shared<psi::PotentialInt>(fact_11.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_22= std::make_shared<psi::PotentialInt>(fact_22.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(1));

  std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
  std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

  potInt_1->set_charge_field(Zxyz_2);
  potInt_2->set_charge_field(Zxyz_1);
  potInt_11->set_charge_field(Zxyz_2);
  potInt_22->set_charge_field(Zxyz_1);
  oneInt = potInt_1;
  oneInt->compute(VaoB12);
  oneInt = potInt_2;
  oneInt->compute(VaoA21);
  oneInt = potInt_11;
  oneInt->compute(VaoB11);
  oneInt = potInt_22;
  oneInt->compute(VaoA22);
  ovlInt->compute(Sao12);
  std::shared_ptr<psi::Matrix> VaoA12 = VaoA21->transpose();
  VaoA21.reset();

  std::shared_ptr<psi::Matrix> Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");

  std::shared_ptr<psi::Matrix> Smo12  = psi::Matrix::triplet(Ca_occ_A, Sao12 , Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA12 = psi::Matrix::triplet(Ca_occ_A, VaoA12, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB12 = psi::Matrix::triplet(Ca_occ_A, VaoB12, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA22 = psi::Matrix::triplet(Ca_occ_B, VaoA22, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB11 = psi::Matrix::triplet(Ca_occ_A, VaoB11, Ca_occ_A, true, false, false);

  VmoA12->add(VmoB12);

  e_s1 = Smo12->vector_dot(VmoA12);
  std::shared_ptr<psi::Matrix> SVS11 = psi::Matrix::triplet(Smo12, VmoA22, Smo12, false,false, true );
  std::shared_ptr<psi::Matrix> SVS22 = psi::Matrix::triplet(Smo12, VmoB11, Smo12, true ,false, false);
  e_s2 = SVS11->trace() + SVS22->trace();

  // ===> Two-electron part <=== //
  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  dpdbuf4 buf_1222, buf_1112, buf_1122;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  double integral;
  double** S = Smo12->pointer();

  global_dpd_->buf4_init(&buf_1112, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,1]"  ), integrals->DPD_ID("[1,2]"  ),
                          integrals->DPD_ID("[1>=1]+"), integrals->DPD_ID("[1,2]"  ), 0, "MO Ints (11|12)");
  global_dpd_->buf4_init(&buf_1222, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[2,2]"  ),
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[2>=2]+"), 0, "MO Ints (12|22)");
  global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,1]"  ), integrals->DPD_ID("[2,2]"  ),
                          integrals->DPD_ID("[1>=1]+"), integrals->DPD_ID("[2>=2]+"), 0, "MO Ints (11|22)");

  int i, j, k, l;
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1112, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1112, h);
       for (int ij = 0; ij < buf_1112.params->rowtot[h]; ++ij) {
            i = buf_1112.params->roworb[h][ij][0];
            j = buf_1112.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1112.params->coltot[h]; ++kl) {
                 k = buf_1112.params->colorb[h][kl][0];
                 l = buf_1112.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 integral = buf_1112.matrix[h][ij][kl]; // ( i_A j_A | k_A l_B )
                 if (i==j) e_s1+= 2.0 * integral * S[k][l];
                 if (j==k) e_s1-=       integral * S[i][l];
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1112, h);
  }
  global_dpd_->buf4_close(&buf_1112);

  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1222, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1222, h);
       for (int ij = 0; ij < buf_1222.params->rowtot[h]; ++ij) {
            i = buf_1222.params->roworb[h][ij][0];
            j = buf_1222.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1222.params->coltot[h]; ++kl) {
                 k = buf_1222.params->colorb[h][kl][0];
                 l = buf_1222.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 integral = buf_1222.matrix[h][ij][kl]; // ( i_A j_B | k_B l_B )
                 if (k==l) e_s1+= 2.0 * integral * S[i][j];
                 if (j==k) e_s1-=       integral * S[i][l];
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1222, h);
  }
  global_dpd_->buf4_close(&buf_1222);

  std::shared_ptr<psi::Matrix> SSmo11 = psi::Matrix::doublet(Smo12, Smo12, false, true );
  std::shared_ptr<psi::Matrix> SSmo22 = psi::Matrix::doublet(Smo12, Smo12, true , false);
  double** SS11 = SSmo11->pointer();
  double** SS22 = SSmo22->pointer();
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
       for (int ij = 0; ij < buf_1122.params->rowtot[h]; ++ij) {
            i = buf_1122.params->roworb[h][ij][0];
            j = buf_1122.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1122.params->coltot[h]; ++kl) {
                 k = buf_1122.params->colorb[h][kl][0];
                 l = buf_1122.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 integral = buf_1122.matrix[h][ij][kl]; // ( i_A j_A | k_B l_B )
                 if (k==l) e_s2+= 2.0 * integral * SS11[i][j];
                 if (i==j) e_s2+= 2.0 * integral * SS22[k][l];
                 e_s2 -= integral * S[i][k] * S[j][l];
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
  }
  global_dpd_->buf4_close(&buf_1122);


  e_s1 *=-2.0;
  e_s2 *= 2.0;

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell et al.)");
  psi::timer_off("Solver E(Paul) Murrell-etal     ");

  // ===> Compute the Exchange Energy <=== //
  double e_exch = compute_pure_exchange_energy();

  // ===> Finish <=== //
  e = e_s1 + e_s2;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (Murrell et al.)             <==\n\n");
     psi::outfile->Printf("     E S^-1    = %13.6f\n", e_s1                             );
     psi::outfile->Printf("     E S^-2    = %13.6f\n", e_s2                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch                           );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch                         );
     psi::outfile->Printf("\n");
  }
 
  // Return the Total Repulsion Energy
  return e; 
}
double RepulsionEnergySolver::compute_benchmark_otto_ladik() {
  double e_s1 = 0.0, e_s2 = 0.0, e = 0;

  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Otto-Ladik)");
  psi::timer_on("Solver E(Paul) Otto-Ladik       ");

  int nbf_1 = wfn_union_->l_nbf(0);
  int nbf_2 = wfn_union_->l_nbf(1);

  // ===> One electron part <=== //
  std::shared_ptr<psi::Matrix> VaoA21    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_1);
  std::shared_ptr<psi::Matrix> VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> VaoA22    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_2);
  std::shared_ptr<psi::Matrix> VaoB11    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)"  , nbf_1, nbf_2);
  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_11(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_22(wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));
  std::shared_ptr<psi::OneBodyAOInt> oneInt, ovlInt(fact_12.ao_overlap());
  std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_12.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_21.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_11= std::make_shared<psi::PotentialInt>(fact_11.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_22= std::make_shared<psi::PotentialInt>(fact_22.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(1));

  std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
  std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

  potInt_1->set_charge_field(Zxyz_2);
  potInt_2->set_charge_field(Zxyz_1);
  potInt_11->set_charge_field(Zxyz_2);
  potInt_22->set_charge_field(Zxyz_1);
  oneInt = potInt_1;
  oneInt->compute(VaoB12);
  oneInt = potInt_2;
  oneInt->compute(VaoA21);
  oneInt = potInt_11;
  oneInt->compute(VaoB11);
  oneInt = potInt_22;
  oneInt->compute(VaoA22);
  ovlInt->compute(Sao12);
  std::shared_ptr<psi::Matrix> VaoA12 = VaoA21->transpose();
  VaoA21.reset();

  std::shared_ptr<psi::Matrix> Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");

  std::shared_ptr<psi::Matrix> Smo12  = psi::Matrix::triplet(Ca_occ_A, Sao12 , Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA12 = psi::Matrix::triplet(Ca_occ_A, VaoA12, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB12 = psi::Matrix::triplet(Ca_occ_A, VaoB12, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoA22 = psi::Matrix::triplet(Ca_occ_B, VaoA22, Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB11 = psi::Matrix::triplet(Ca_occ_A, VaoB11, Ca_occ_A, true, false, false);

  VmoA12->add(VmoB12);

  e_s1 = Smo12->vector_dot(VmoA12);

  double** S = Smo12->pointer();
  for (int a=0; a<wfn_union_->l_ndocc(0); ++a) {
       for (int b=0; b<wfn_union_->l_ndocc(1); ++b) {
            e_s2 += S[a][b] * S[a][b] * (VmoB11->get(a,a) + VmoA22->get(b,b));
       }
  }

  // ===> Two-electron part <=== //
  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  dpdbuf4 buf_1222, buf_1112, buf_1122;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  double integral;

  global_dpd_->buf4_init(&buf_1112, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,1]"  ), integrals->DPD_ID("[1,2]"  ),
                          integrals->DPD_ID("[1>=1]+"), integrals->DPD_ID("[1,2]"  ), 0, "MO Ints (11|12)");
  global_dpd_->buf4_init(&buf_1222, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[2,2]"  ),
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[2>=2]+"), 0, "MO Ints (12|22)");
  global_dpd_->buf4_init(&buf_1122, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,1]"  ), integrals->DPD_ID("[2,2]"  ),
                          integrals->DPD_ID("[1>=1]+"), integrals->DPD_ID("[2>=2]+"), 0, "MO Ints (11|22)");

  int i, j, k, l;
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1112, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1112, h);
       for (int ij = 0; ij < buf_1112.params->rowtot[h]; ++ij) {
            i = buf_1112.params->roworb[h][ij][0];
            j = buf_1112.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1112.params->coltot[h]; ++kl) {
                 k = buf_1112.params->colorb[h][kl][0];
                 l = buf_1112.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 if (i==j) {
                           integral = buf_1112.matrix[h][ij][kl]; // ( i_A j_A | k_A l_B )
                           e_s1+= 2.0 * integral * S[k][l];
                 if (i==k) e_s1-=       integral * S[k][l];
                 }
                 
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1112, h);
  }
  global_dpd_->buf4_close(&buf_1112);

  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1222, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1222, h);
       for (int ij = 0; ij < buf_1222.params->rowtot[h]; ++ij) {
            i = buf_1222.params->roworb[h][ij][0];
            j = buf_1222.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1222.params->coltot[h]; ++kl) {
                 k = buf_1222.params->colorb[h][kl][0];
                 l = buf_1222.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 if (k==l) {
                           integral = buf_1222.matrix[h][ij][kl]; // ( i_A j_B | k_B l_B )
                           e_s1+= 2.0 * integral * S[i][j];
                 if (j==k) e_s1-=       integral * S[i][j];
                 }
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1222, h);
  }
  global_dpd_->buf4_close(&buf_1222);

  std::shared_ptr<psi::Vector> ss_1 = std::make_shared<psi::Vector>("s^2 (1)", wfn_union_->l_ndocc(0));
  std::shared_ptr<psi::Vector> ss_2 = std::make_shared<psi::Vector>("s^2 (2)", wfn_union_->l_ndocc(1));
  double* s1 = ss_1->pointer();
  double* s2 = ss_2->pointer();
  double val;
  for (int a = 0; a<wfn_union_->l_ndocc(0); ++a) {
       val = 0.0;
       for (int b = 0; b<wfn_union_->l_ndocc(1); ++b) {
            val += S[a][b] * S[a][b];
       }
       s1[a] = val;
  }
  for (int b = 0; b<wfn_union_->l_ndocc(1); ++b) {
       val = 0.0;
       for (int a = 0; a<wfn_union_->l_ndocc(0); ++a) {
            val += S[a][b] * S[a][b];
       }
       s2[b] = val;
  }

  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1122, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1122, h);
       for (int ij = 0; ij < buf_1122.params->rowtot[h]; ++ij) {
            i = buf_1122.params->roworb[h][ij][0];
            j = buf_1122.params->roworb[h][ij][1];
            for (int kl = 0; kl < buf_1122.params->coltot[h]; ++kl) {
                 k = buf_1122.params->colorb[h][kl][0];
                 l = buf_1122.params->colorb[h][kl][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 if ((i==j) && (k==l)) {
                     integral = buf_1122.matrix[h][ij][kl]; // ( i_A j_A | k_B l_B )
                     e_s2 += 2.0 * integral * (s1[i] + s2[k]); 
                     e_s2 -=       integral * S[i][k] * S[i][k];
                 }
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1122, h);
  }
  global_dpd_->buf4_close(&buf_1122);


  e_s1 *=-2.0;
  e_s2 *= 2.0;

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Otto-Ladik)");
  psi::timer_off("Solver E(Paul) Otto-Ladik       ");

  // ===> Compute the Exchange Energy <=== //
  double e_exch = compute_pure_exchange_energy();

  // ===> Finish <=== //
  e = e_s1 + e_s2;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (Otto-Ladik)                 <==\n\n");
     psi::outfile->Printf("     E S^-1    = %13.6f\n", e_s1                             );
     psi::outfile->Printf("     E S^-2    = %13.6f\n", e_s2                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch                           );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch                         );
     psi::outfile->Printf("\n");
  }
 
  // Return the Total Repulsion Energy
  return e; 
}
double RepulsionEnergySolver::compute_benchmark_efp2() {
  double e_s1 = 0.0, e_s2 = 0.0, e = 0;

  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (EFP2)");
  psi::timer_on ("Solver E(Paul) EFP2             ");

  int nbf_1 = wfn_union_->l_nbf(0);
  int nbf_2 = wfn_union_->l_nbf(1);
  std::shared_ptr<psi::Molecule> mol_1 = wfn_union_->l_wfn(0)->molecule();
  std::shared_ptr<psi::Molecule> mol_2 = wfn_union_->l_wfn(1)->molecule();
  int nat_1 = mol_1->natom();
  int nat_2 = mol_2->natom();

  // Compute orbital centroids as well as overlap, kinetic and Fock matrix elements in MO basis
  std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)", nbf_1, nbf_2);
  std::shared_ptr<psi::Matrix> Tao12     = std::make_shared<psi::Matrix>("Tao(1,2)", nbf_1, nbf_2);
  std::vector<std::shared_ptr<psi::Matrix>> R1ao, R2ao;
  std::vector<std::shared_ptr<psi::Vector>> R1mo, R2mo;

  for (int z=0; z<3; ++z) R1ao.push_back(std::make_shared<psi::Matrix>("R1ao", nbf_1, nbf_1));
  for (int z=0; z<3; ++z) R2ao.push_back(std::make_shared<psi::Matrix>("R2ao", nbf_2, nbf_2));

  for (int z=0; z<3; ++z) R1mo.push_back(std::make_shared<psi::Vector>("R1mo", wfn_union_->l_ndocc(0)));
  for (int z=0; z<3; ++z) R2mo.push_back(std::make_shared<psi::Vector>("R2mo", wfn_union_->l_ndocc(1)));

  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1), wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_11(wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_22(wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1), wfn_union_->l_primary(1));

  std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_12.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> kinInt(fact_12.ao_kinetic());
  std::shared_ptr<psi::OneBodyAOInt> dipInt1(fact_11.ao_dipole());
  std::shared_ptr<psi::OneBodyAOInt> dipInt2(fact_22.ao_dipole());

  ovlInt->compute(Sao12);
  kinInt->compute(Tao12);
  dipInt1->compute(R1ao);
  dipInt2->compute(R2ao);

  std::shared_ptr<psi::Matrix> Ca_occ_A = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_B = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");
  for (int z=0; z<3; ++z) {
       R1ao[z]->scale(-1.0);
       R2ao[z]->scale(-1.0);
       std::shared_ptr<psi::Matrix> CR1C = psi::Matrix::triplet(Ca_occ_A, R1ao[z], Ca_occ_A, true, false, false);
       std::shared_ptr<psi::Matrix> CR2C = psi::Matrix::triplet(Ca_occ_B, R2ao[z], Ca_occ_B, true, false, false);
       for (int a=0; a<wfn_union_->l_ndocc(0); ++a) {
            R1mo[z]->set(a, CR1C->get(a,a));
       }
       for (int b=0; b<wfn_union_->l_ndocc(1); ++b) {
            R2mo[z]->set(b, CR2C->get(b,b));
       }
       R1ao[z].reset(); 
       R2ao[z].reset(); 
  }

  std::shared_ptr<psi::Matrix> Smo12 = psi::Matrix::triplet(Ca_occ_A, Sao12 , Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> Tmo12 = psi::Matrix::triplet(Ca_occ_A, Tao12 , Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> F1    = psi::Matrix::triplet(Ca_occ_A, wfn_union_->l_wfn(0)->Fa(), Ca_occ_A, true, false, false);
  std::shared_ptr<psi::Matrix> F2    = psi::Matrix::triplet(Ca_occ_B, wfn_union_->l_wfn(1)->Fa(), Ca_occ_B, true, false, false);
  std::shared_ptr<psi::Matrix> SF1S  = psi::Matrix::triplet(Smo12, F1, Smo12, true, false, false);
  std::shared_ptr<psi::Matrix> SF2S  = psi::Matrix::triplet(Smo12, F2, Smo12, false, false, true);

  e_s1 = SF1S->trace() + SF2S->trace() - 2.0 * Smo12->vector_dot(Tmo12);

  SF1S.reset();
  SF2S.reset();

  double** S = Smo12->pointer();
  std::shared_ptr<psi::Vector> ss1 = std::make_shared<psi::Vector>("s^2 (1)", wfn_union_->l_ndocc(0));
  std::shared_ptr<psi::Vector> ss2 = std::make_shared<psi::Vector>("s^2 (2)", wfn_union_->l_ndocc(1));
  double* s1 = ss1->pointer();
  double* s2 = ss2->pointer();
  double val;
  for (int a = 0; a<wfn_union_->l_ndocc(0); ++a) {
       val = 0.0;
       for (int b = 0; b<wfn_union_->l_ndocc(1); ++b) {
            val += S[a][b] * S[a][b];
       }
       s1[a] = val;
  }
  for (int b = 0; b<wfn_union_->l_ndocc(1); ++b) {
       val = 0.0;
       for (int a = 0; a<wfn_union_->l_ndocc(0); ++a) {
            val += S[a][b] * S[a][b];
       }
       s2[b] = val;
  }

  std::shared_ptr<psi::Vector> ZR1 = std::make_shared<psi::Vector>("ZR1", wfn_union_->l_ndocc(0));
  std::shared_ptr<psi::Vector> ZR2 = std::make_shared<psi::Vector>("ZR2", wfn_union_->l_ndocc(1));
  std::shared_ptr<psi::Vector> IR1 = std::make_shared<psi::Vector>("IR1", wfn_union_->l_ndocc(0));
  std::shared_ptr<psi::Vector> IR2 = std::make_shared<psi::Vector>("IR2", wfn_union_->l_ndocc(1));

  double* pZR1 = ZR1->pointer();
  double* pZR2 = ZR2->pointer();
  double* pIR1 = IR1->pointer();
  double* pIR2 = IR2->pointer();
 
  double rbx, ray, rad, rbc, rab;

  for (int b=0; b<wfn_union_->l_ndocc(1); ++b){
       val = 0.0;
       for (int x=0; x<nat_1; ++x) { 
            rbx  = sqrt(pow(mol_1->x(x)-R2mo[0]->get(b), 2.0) + 
                        pow(mol_1->y(x)-R2mo[1]->get(b), 2.0) + 
                        pow(mol_1->z(x)-R2mo[2]->get(b), 2.0) );
            val -= (double)mol_1->Z(x) / rbx;
       }
       pZR2[b] = val;
  }
  for (int a=0; a<wfn_union_->l_ndocc(0); ++a){
       val = 0.0;
       for (int y=0; y<nat_2; ++y) { 
            ray  = sqrt(pow(mol_2->x(y)-R1mo[0]->get(a), 2.0) + 
                        pow(mol_2->y(y)-R1mo[1]->get(a), 2.0) + 
                        pow(mol_2->z(y)-R1mo[2]->get(a), 2.0) );
            val -= (double)mol_2->Z(y) / ray;
       }
       pZR1[a] = val;
  }
  for (int a=0; a<wfn_union_->l_ndocc(0); ++a) {
       val = 0.0;
       for (int b=0; b<wfn_union_->l_ndocc(1); ++b) {
            rad  = sqrt(pow(R1mo[0]->get(a)-R2mo[0]->get(b), 2.0) + 
                        pow(R1mo[1]->get(a)-R2mo[1]->get(b), 2.0) + 
                        pow(R1mo[2]->get(a)-R2mo[2]->get(b), 2.0) );
            val += 2.0 / rad;
       }
       pIR1[a] = val;
  }
  for (int b=0; b<wfn_union_->l_ndocc(1); ++b) {
       val = 0.0;
       for (int a=0; a<wfn_union_->l_ndocc(0); ++a) {
            rbc  = sqrt(pow(R1mo[0]->get(a)-R2mo[0]->get(b), 2.0) + 
                        pow(R1mo[1]->get(a)-R2mo[1]->get(b), 2.0) + 
                        pow(R1mo[2]->get(a)-R2mo[2]->get(b), 2.0) );
            val += 2.0 / rbc;
       }
       pIR2[b] = val;
  }
  ZR1->add(IR1);
  ZR2->add(IR2);

  e_s2 = ss1->vector_dot(ZR1) + ss2->vector_dot(ZR2);

  for (int a=0; a<wfn_union_->l_ndocc(0); ++a) {
       for (int b=0; b<wfn_union_->l_ndocc(1); ++b) { 
            rab = sqrt(pow(R1mo[0]->get(a) - R2mo[0]->get(b), 2.0) +
                       pow(R1mo[1]->get(a) - R2mo[1]->get(b), 2.0) +
                       pow(R1mo[2]->get(a) - R2mo[2]->get(b), 2.0) );
            e_s2 -= S[a][b] * S[a][b] / rab;
       }
  }

  // ---> Finalize with repulsion <--- //
  e_s1 *=-2.0;
  e_s2 *= 2.0;
  //psi::timer_off("SOLVER: Repulsion Energy Calculations (EFP2)");
  psi::timer_off("Solver E(Paul) EFP2             ");

  // ===> Compute the Exchange Energy <=== //
  double e_exch_pure = compute_pure_exchange_energy();
  double e_exch_efp2 = compute_efp2_exchange_energy(Smo12, R1mo, R2mo);

  // ===> Finish <=== //
  e = e_s1 + e_s2;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>         Benchmark (EFP2)                       <==\n\n");
     psi::outfile->Printf("     E S^-1    = %13.6f\n", e_s1                             );
     psi::outfile->Printf("     E S^-2    = %13.6f\n", e_s2                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch_pure                      );
     psi::outfile->Printf("     E EX(SGO) = %13.6f\n", e_exch_efp2                      );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch_pure                    );
     psi::outfile->Printf("     EXREP(SGO)= %13.6f\n", e+e_exch_efp2                    );
     psi::outfile->Printf("\n");
  }
 
  // Return the Total Repulsion Energy
  return e; 
}
double RepulsionEnergySolver::compute_efp2_exchange_energy(psi::SharedMatrix S, 
          std::vector<psi::SharedVector> R1, std::vector<psi::SharedVector> R2) {

 //psi::timer_on ("SOLVER: Exchange Energy Calculations (SGO)");
 psi::timer_on("Solver E(Exch) EFP2(SGO)        ");

 double e = 0.0, sab, rab;

 double** s = S->pointer();
 for (int a=0; a<S->nrow(); ++a) {
      for (int b=0; b<S->ncol(); ++b) {
           sab = s[a][b];
           rab = sqrt(pow(R1[0]->get(a) - R2[0]->get(b), 2.0) + 
                      pow(R1[1]->get(a) - R2[1]->get(b), 2.0) + 
                      pow(R1[2]->get(a) - R2[2]->get(b), 2.0) );
           e += sqrt(-2.0*log(abs(sab)) / M_PI) * sab * sab / rab;
      }
 }
 e *= -4.0;
 //psi::timer_off("SOLVER: Exchange Energy Calculations (SGO)");
 psi::timer_off("Solver E(Exch) EFP2(SGO)        ");

 return e;
}
double RepulsionEnergySolver::compute_oep_based_murrell_etal_gdf_camm() {
  double e = 0.0, e_s1 = 0.0, e_s2 = 0.0;

  SharedOEPotential oep_1 = oepdev::OEPotential::build("REPULSION ENERGY",
                                                       wfn_union_->l_wfn(0), 
                                                       wfn_union_->l_auxiliary(0), 
                                                       wfn_union_->l_intermediate(0), 
                                                       wfn_union_->options());
  SharedOEPotential oep_2 = oepdev::OEPotential::build("REPULSION ENERGY", 
                                                       wfn_union_->l_wfn(1), 
                                                       wfn_union_->l_auxiliary(1), 
                                                       wfn_union_->l_intermediate(1), 
                                                       wfn_union_->options());
  oep_1->compute();
  oep_2->compute();

  // ===> Compute Overlap Matrices <=== //
  int nbf_p1 = wfn_union_->l_nbf(0);
  int nbf_p2 = wfn_union_->l_nbf(1);
  int nbf_a1 = wfn_union_->l_auxiliary(0)->nbf();
  int nbf_a2 = wfn_union_->l_auxiliary(1)->nbf();

  std::shared_ptr<psi::Matrix> Sao_1p2p     = std::make_shared<psi::Matrix>("Sao 1p2p", nbf_p1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1a2p     = std::make_shared<psi::Matrix>("Sao 1a2p", nbf_a1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1p2a     = std::make_shared<psi::Matrix>("Sao 1p2a", nbf_p1, nbf_a2);

  psi::IntegralFactory fact_1p2p(wfn_union_->l_primary  (0), wfn_union_->l_primary  (1), wfn_union_->l_primary  (0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1a2p(wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1), wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1p2a(wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1), wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1));

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2p(fact_1p2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

  ovlInt_1p2p->compute(Sao_1p2p);
  ovlInt_1a2p->compute(Sao_1a2p);
  ovlInt_1p2a->compute(Sao_1p2a);

  std::shared_ptr<psi::Matrix> Ca_occ_1 = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_2 = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");

  std::shared_ptr<psi::Matrix> Smo = psi::Matrix::triplet(Ca_occ_1, Sao_1p2p, Ca_occ_2, true, false, false);
  std::shared_ptr<psi::Matrix> Sba = psi::Matrix::doublet(Ca_occ_2, Sao_1a2p, true, true);
  std::shared_ptr<psi::Matrix> Sab = psi::Matrix::doublet(Ca_occ_1, Sao_1p2a, true, false);

  // ===> Compute S^-1 term <=== //
  std::shared_ptr<psi::Matrix> SSG1= psi::Matrix::triplet(Smo, Sba, oep_1->matrix("Murrell-etal.S1"), 
                                                          false, false, false);
  std::shared_ptr<psi::Matrix> SSG2= psi::Matrix::triplet(Smo, Sab, oep_2->matrix("Murrell-etal.S1"), 
                                                          true, false, false);
  e_s1  = SSG1->trace() + SSG2->trace();
  e_s1 *= -2.0;

  // ===> Compute S^-2 term <=== //
  SharedMTPConv conv_aB = oep_1->oep("Otto-Ladik.S2.CAMM.a").dmtp->energy(oep_2->oep("Otto-Ladik.S2.CAMM.A").dmtp);
  SharedMTPConv conv_bA = oep_2->oep("Otto-Ladik.S2.CAMM.a").dmtp->energy(oep_1->oep("Otto-Ladik.S2.CAMM.A").dmtp);

  MultipoleConvergence::ConvergenceLevel clevel;
  if      (options_.get_str("DMTP_CONVER") == "R1") clevel = MultipoleConvergence::ConvergenceLevel::R1;
  else if (options_.get_str("DMTP_CONVER") == "R2") clevel = MultipoleConvergence::ConvergenceLevel::R2;
  else if (options_.get_str("DMTP_CONVER") == "R3") clevel = MultipoleConvergence::ConvergenceLevel::R3;
  else if (options_.get_str("DMTP_CONVER") == "R4") clevel = MultipoleConvergence::ConvergenceLevel::R4;
  else if (options_.get_str("DMTP_CONVER") == "R5") clevel = MultipoleConvergence::ConvergenceLevel::R5;
  else {throw psi::PSIEXCEPTION("Incorrect convergence level specified!");}

  psi::SharedMatrix level_aB = conv_aB->level(clevel);
  psi::SharedMatrix level_bA = conv_bA->level(clevel);

  for (int i=0; i< oep_1->n("Otto-Ladik.S2.CAMM.a"); ++i) {
  for (int j=0; j< oep_2->n("Otto-Ladik.S2.CAMM.a"); ++j) {
      double sij = Smo->get(i,j);
      e_s2 += sij * sij * (level_aB->get(i,j) + level_bA->get(j,i));
  }
  }
  e_s2 *= 2.0;
  
  // ===> Compute the Exchange Energy <=== //
  double e_exch_pure = compute_pure_exchange_energy();
  double e_exch_efp2 = 0.0;//compute_efp2_exchange_energy(Smo12, R1mo, R2mo);

  // ===> Finish <=== //
  e = e_s1 + e_s2;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>    OEP-Based (Murrell-etal; S1-GDF/S2-CAMM)    <==\n\n");
     psi::outfile->Printf("     E S^-1    = %13.6f\n", e_s1                             );
     psi::outfile->Printf("     E S^-2    = %13.6f\n", e_s2                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch_pure                      );
   //psi::outfile->Printf("     E EX(SGO) = %13.6f\n", e_exch_efp2                      );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch_pure                    );
   //psi::outfile->Printf("     EXREP(SGO)= %13.6f\n", e+e_exch_efp2                    );
     psi::outfile->Printf("\n");
  }
 
  // Return the Total Repulsion Energy
  return e; 

}
double RepulsionEnergySolver::compute_oep_based_murrell_etal_gdf_esp() {
  double e = 0.0, e_s1 = 0.0, e_s2 = 0.0;

  SharedOEPotential oep_1 = oepdev::OEPotential::build("REPULSION ENERGY",
                                                       wfn_union_->l_wfn(0), 
                                                       wfn_union_->l_auxiliary(0), 
                                                       wfn_union_->l_intermediate(0), 
                                                       wfn_union_->options());
  SharedOEPotential oep_2 = oepdev::OEPotential::build("REPULSION ENERGY", 
                                                       wfn_union_->l_wfn(1), 
                                                       wfn_union_->l_auxiliary(1), 
                                                       wfn_union_->l_intermediate(1), 
                                                       wfn_union_->options());
  oep_1->compute();
  oep_2->compute();

  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell-OEP:S1-DF/S2-ESP)");
  psi::timer_on("Solver E(Paul) OEP-Based        ");

  // ===> Compute Overlap Matrices <=== //
  int nbf_p1 = wfn_union_->l_nbf(0);
  int nbf_p2 = wfn_union_->l_nbf(1);
  int nbf_a1 = wfn_union_->l_auxiliary(0)->nbf();
  int nbf_a2 = wfn_union_->l_auxiliary(1)->nbf();

  std::shared_ptr<psi::Matrix> Sao_1p2p     = std::make_shared<psi::Matrix>("Sao 1p2p", nbf_p1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1a2p     = std::make_shared<psi::Matrix>("Sao 1a2p", nbf_a1, nbf_p2);
  std::shared_ptr<psi::Matrix> Sao_1p2a     = std::make_shared<psi::Matrix>("Sao 1p2a", nbf_p1, nbf_a2);

  psi::IntegralFactory fact_1p2p(wfn_union_->l_primary  (0), wfn_union_->l_primary  (1), wfn_union_->l_primary  (0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1a2p(wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1), wfn_union_->l_auxiliary(0), wfn_union_->l_primary  (1));
  psi::IntegralFactory fact_1p2a(wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1), wfn_union_->l_primary  (0), wfn_union_->l_auxiliary(1));

  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2p(fact_1p2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1a2p(fact_1a2p.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> ovlInt_1p2a(fact_1p2a.ao_overlap());

  ovlInt_1p2p->compute(Sao_1p2p);
  ovlInt_1a2p->compute(Sao_1a2p);
  ovlInt_1p2a->compute(Sao_1p2a);

  std::shared_ptr<psi::Matrix> Ca_occ_1 = wfn_union_->l_wfn(0)->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_occ_2 = wfn_union_->l_wfn(1)->Ca_subset("AO","OCC");

  std::shared_ptr<psi::Matrix> Smo = psi::Matrix::triplet(Ca_occ_1, Sao_1p2p, Ca_occ_2, true, false, false);
  std::shared_ptr<psi::Matrix> Sba = psi::Matrix::doublet(Ca_occ_2, Sao_1a2p, true, true);
  std::shared_ptr<psi::Matrix> Sab = psi::Matrix::doublet(Ca_occ_1, Sao_1p2a, true, false);

  // ===> Compute S^-1 term <=== //
  std::shared_ptr<psi::Matrix> SSG1= psi::Matrix::triplet(Smo, Sba, oep_1->matrix("Murrell-etal.S1"), 
                                                          false, false, false);
  std::shared_ptr<psi::Matrix> SSG2= psi::Matrix::triplet(Smo, Sab, oep_2->matrix("Murrell-etal.S1"), 
                                                          true, false, false);
  e_s1  = SSG1->trace() + SSG2->trace();

  // ===> Compute S^-2 term <=== //
  std::vector<std::shared_ptr<psi::Matrix>> V1_set;
  std::vector<std::shared_ptr<psi::Matrix>> V2_set;
  psi::IntegralFactory fact_1p1p(wfn_union_->l_primary  (0), wfn_union_->l_primary(0), wfn_union_->l_primary  (0), wfn_union_->l_primary(0));
  psi::IntegralFactory fact_2p2p(wfn_union_->l_primary  (1), wfn_union_->l_primary(1), wfn_union_->l_primary  (1), wfn_union_->l_primary(1));

  for (int nx=0; nx<oep_1->oep("Otto-Ladik.S2.ESP").n; ++nx) {
       std::shared_ptr<psi::Matrix> V2 = std::make_shared<psi::Matrix>("V2", nbf_p2, nbf_p2);

       std::shared_ptr<psi::OneBodyAOInt> oneInt;
       std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_2p2p.spherical_transform(),
                                                                                         wfn_union_->l_primary(1),
                                                                                         wfn_union_->l_primary(1));

       std::shared_ptr<psi::Matrix> Qxyz_1 = std::make_shared<psi::Matrix>("Q OEP 1", oep_1->matrix("Otto-Ladik.S2.ESP")->nrow(), 4); 
       for (int i=0; i < Qxyz_1->nrow(); ++i) {
            Qxyz_1->set(i, 0, oep_1->matrix("Otto-Ladik.S2.ESP")->get(i,nx));
            Qxyz_1->set(i, 1, oep_1->wfn()->molecule()->x(i));
            Qxyz_1->set(i, 2, oep_1->wfn()->molecule()->y(i));
            Qxyz_1->set(i, 3, oep_1->wfn()->molecule()->z(i));
       }
       potInt_2->set_charge_field(Qxyz_1);
       oneInt = potInt_2;
       oneInt->compute(V2);
       //
       V2_set.push_back(psi::Matrix::triplet(Ca_occ_2, V2, Ca_occ_2, true, false, false));
  }

  for (int ny=0; ny<oep_2->oep("Otto-Ladik.S2.ESP").n; ++ny) {
       std::shared_ptr<psi::Matrix> V1 = std::make_shared<psi::Matrix>("V1", nbf_p1, nbf_p1);

       std::shared_ptr<psi::OneBodyAOInt> oneInt;
       std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_1p1p.spherical_transform(),
                                                                                         wfn_union_->l_primary(0),
                                                                                         wfn_union_->l_primary(0));

       std::shared_ptr<psi::Matrix> Qxyz_2 = std::make_shared<psi::Matrix>("Q OEP 2", oep_2->matrix("Otto-Ladik.S2.ESP")->nrow(), 4); 
       for (int i=0; i < Qxyz_2->nrow(); ++i) {
            Qxyz_2->set(i, 0, oep_2->matrix("Otto-Ladik.S2.ESP")->get(i,ny));
            Qxyz_2->set(i, 1, oep_2->wfn()->molecule()->x(i));
            Qxyz_2->set(i, 2, oep_2->wfn()->molecule()->y(i));
            Qxyz_2->set(i, 3, oep_2->wfn()->molecule()->z(i));
       }
       potInt_1->set_charge_field(Qxyz_2);
       oneInt = potInt_1;
       oneInt->compute(V1);
       //
       V1_set.push_back(psi::Matrix::triplet(Ca_occ_1, V1, Ca_occ_1, true, false, false));
  }

  for (int ox=0; ox< oep_1->oep("Otto-Ladik.S2.ESP").n; ++ox) {
  for (int oy=0; oy< oep_2->oep("Otto-Ladik.S2.ESP").n; ++oy) {
       double v = Smo->get(ox, oy);
       e_s2 += v * v * (V2_set[ox]->get(oy, oy) + V1_set[oy]->get(ox, ox));
  }
  }

  e_s1 *= -2.0;
  e_s2 *=  2.0; 

  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell-OEP:S1-DF/S2-ESP)");
  psi::timer_off("Solver E(Paul) OEP-Based        ");

  // ===> Compute the Exchange Energy <=== //
  double e_exch_pure = compute_pure_exchange_energy();
  double e_exch_efp2 = 0.0;//compute_efp2_exchange_energy(Smo12, R1mo, R2mo);

  // ===> Finish <=== //
  e = e_s1 + e_s2;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n"  );
     psi::outfile->Printf("  ==>     OEP-Based (Murrell-etal; S1-GDF/S2-ESP)    <==\n\n");
     psi::outfile->Printf("     E S^-1    = %13.6f\n", e_s1                             );
     psi::outfile->Printf("     E S^-2    = %13.6f\n", e_s2                             );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     E REP     = %13.6f\n", e                                );
     psi::outfile->Printf("     E EX      = %13.6f\n", e_exch_pure                      );
   //psi::outfile->Printf("     E EX(SGO) = %13.6f\n", e_exch_efp2                      );
     psi::outfile->Printf("     -------------------------------\n"                      );
     psi::outfile->Printf("     EXREP     = %13.6f\n", e+e_exch_pure                    );
   //psi::outfile->Printf("     EXREP(SGO)= %13.6f\n", e+e_exch_efp2                    );
     psi::outfile->Printf("\n");
  }
 
  // Return the Total Repulsion Energy
  return e; 
}
double RepulsionEnergySolver::compute_oep_based_murrell_etal_esp() {
  //psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell-OEP:ESP)");
  //psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell-OEP:ESP)");
  throw psi::PSIEXCEPTION("ERROR: MURRELL_ETAL_ESP is not yet implemented!\n");
}
// Build: factory static method 
std::shared_ptr<OEPDevSolver> OEPDevSolver::build(const std::string& target, SharedWavefunctionUnion wfn_union)
{
   std::shared_ptr<OEPDevSolver> solver;
   if      (target == "ELECTROSTATIC ENERGY"  ) solver = std::make_shared< ElectrostaticEnergySolver>(wfn_union);
   else if (target == "REPULSION ENERGY"      ) solver = std::make_shared<     RepulsionEnergySolver>(wfn_union);
   else if (target == "CHARGE TRANSFER ENERGY") solver = std::make_shared<ChargeTransferEnergySolver>(wfn_union);
   else throw psi::PSIEXCEPTION("OEPDEV: Error. OEPDevSolver: Incorrect targer property chosen!\n");
   return solver;
}
