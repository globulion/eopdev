//#include "psi4/libtrans/integraltransform.h"
//#include "psi4/libdpd/dpd.h"

#include "solver.h"

using namespace std;
using namespace psi;
using namespace oepdev;

OEPDevSolver::OEPDevSolver(SharedWavefunctionUnion wfn_union)
 : wfn_union_(wfn_union), methods_oepBased_(), methods_benchmark_()
{
  
}
OEPDevSolver::~OEPDevSolver()
{

}
double OEPDevSolver::compute_oep_based(const std::string& method) {}
double OEPDevSolver::compute_benchmark(const std::string& method) {}
ElectrostaticEnergySolver::ElectrostaticEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  methods_oepBased_ .push_back("ESP_SYMMETRIZED");
  methods_benchmark_.push_back("AO_EXPANDED"    );
  methods_benchmark_.push_back("MO_EXPANDED"    );
}
ElectrostaticEnergySolver::~ElectrostaticEnergySolver() 
{

}
double ElectrostaticEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  double e_nuc_mol = 0.0;
  double e_mol_nuc = 0.0;
  double e_el_mol  = 0.0;
  double e_mol_el  = 0.0;
  double qq        = 0.0;

  if (method == "DEFAULT" || method == "ESP_SYMMETRIZED") {
 
     // ===> [Nuc+El](A) --- El(B) <=== //
     SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
     SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
     SharedOEPotential oep_1 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_1, wfn_union_->options());
     SharedOEPotential oep_2 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_2, wfn_union_->options());

     oep_1->compute();
     oep_2->compute();

     oep_1->print_header();
     oep_2->print_header();

     psi::timer_on("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");

     int nbf_1 = wfn_union_->l_primary(0)->nbf();
     int nbf_2 = wfn_union_->l_primary(1)->nbf();

     std::shared_ptr<psi::Matrix> V1 = std::make_shared<psi::Matrix>("V1", nbf_1, nbf_1);
     std::shared_ptr<psi::Matrix> V2 = std::make_shared<psi::Matrix>("V2", nbf_2, nbf_2);

     psi::IntegralFactory fact_1(wfn_union_->l_primary(0), wfn_union_->l_primary(0));
     psi::IntegralFactory fact_2(wfn_union_->l_primary(1), wfn_union_->l_primary(1));

     std::shared_ptr<psi::OneBodyAOInt> oneInt;
     std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_1.spherical_transform(),
                                                                                       wfn_union_->l_primary(0),
                                                                                       wfn_union_->l_primary(0));
     std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_2.spherical_transform(),
                                                                                       wfn_union_->l_primary(1),
                                                                                       wfn_union_->l_primary(1));

     std::shared_ptr<psi::Matrix> Qxyz_1 = std::make_shared<psi::Matrix>("Q OEP 1", oep_1->matrix("V")->nrow(), 4);
     std::shared_ptr<psi::Matrix> Qxyz_2 = std::make_shared<psi::Matrix>("Q OEP 2", oep_2->matrix("V")->nrow(), 4);
     for (int i=0; i < Qxyz_1->nrow(); ++i) {
          Qxyz_1->set(i, 0, oep_1->matrix("V")->get(i,0));
          Qxyz_1->set(i, 1, oep_1->wfn()->molecule()->x(i));
          Qxyz_1->set(i, 2, oep_1->wfn()->molecule()->y(i));
          Qxyz_1->set(i, 3, oep_1->wfn()->molecule()->z(i));
     }
     for (int i=0; i < Qxyz_2->nrow(); ++i) {
          Qxyz_2->set(i, 0, oep_2->matrix("V")->get(i,0));
          Qxyz_2->set(i, 1, oep_2->wfn()->molecule()->x(i));
          Qxyz_2->set(i, 2, oep_2->wfn()->molecule()->y(i));
          Qxyz_2->set(i, 3, oep_2->wfn()->molecule()->z(i));
     }

     // ===> ESP(A) --- NUC(B) <=== //
     for (int i=0; i < Qxyz_1->nrow(); ++i) {
     for (int j=0; j < oep_2->wfn()->molecule()->natom(); ++j) {
          e_mol_nuc += Qxyz_1->get(i, 0) * oep_2->wfn()->molecule()->Z(j) /
               sqrt( pow( Qxyz_1->get(i,1) - oep_2->wfn()->molecule()->x(j), 2.0) +
                     pow( Qxyz_1->get(i,2) - oep_2->wfn()->molecule()->y(j), 2.0) +
                     pow( Qxyz_1->get(i,3) - oep_2->wfn()->molecule()->z(j), 2.0) );
     }
     }
     // ===> Nuc(A) --- ESP(B) <=== //
     double qe = 0.0;
     for (int i=0; i < oep_1->wfn()->molecule()->natom(); ++i) {
     for (int j=0; j < Qxyz_2->nrow(); ++j) {
          e_nuc_mol += Qxyz_2->get(j, 0) * oep_1->wfn()->molecule()->Z(i) /
               sqrt( pow( Qxyz_2->get(j,1) - oep_1->wfn()->molecule()->x(i), 2.0) +
                     pow( Qxyz_2->get(j,2) - oep_1->wfn()->molecule()->y(i), 2.0) +
                     pow( Qxyz_2->get(j,3) - oep_1->wfn()->molecule()->z(i), 2.0) );
     }
     }

     potInt_1->set_charge_field(Qxyz_2);
     potInt_2->set_charge_field(Qxyz_1);

     oneInt = potInt_1;
     oneInt->compute(V1);
     oneInt = potInt_2;
     oneInt->compute(V2);

     e_el_mol += wfn_1->Da()->vector_dot(V1);
     e_el_mol += wfn_1->Db()->vector_dot(V1);
     e_mol_el += wfn_2->Da()->vector_dot(V2);
     e_mol_el += wfn_2->Db()->vector_dot(V2);

     // Finish
     e = 0.5*(e_mol_nuc + e_mol_el + e_nuc_mol + e_el_mol);

     psi::timer_off("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");

     // ===> ESP(A) ---- ESP(B) <=== //
     psi::timer_on ("SOLVER: Electrostatic Energy Calculations (ESP-ESP)");

     for (int i=0; i < Qxyz_1->nrow(); ++i) {
     for (int j=0; j < Qxyz_2->nrow(); ++j) {
          qq += Qxyz_1->get(i, 0) * Qxyz_2->get(j, 0) /
               sqrt( pow( Qxyz_1->get(i,1) - Qxyz_2->get(j,1), 2.0) +
                     pow( Qxyz_1->get(i,2) - Qxyz_2->get(j,2), 2.0) +
                     pow( Qxyz_1->get(i,3) - Qxyz_2->get(j,3), 2.0) );
     }
     }
     psi::timer_off("SOLVER: Electrostatic Energy Calculations (ESP-ESP)");


     // Print
     if (wfn_union_->options().get_int("PRINT") > 0) {
        psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
        psi::outfile->Printf("  ==>         OEP-based Model (DEFAULT)         <==\n\n");
        psi::outfile->Printf("     ESP NUC   = %13.6f\n", e_mol_nuc);
        psi::outfile->Printf("     ESP EL    = %13.6f\n", e_mol_el );
        psi::outfile->Printf("     NUC ESP   = %13.6f\n", e_nuc_mol);
        psi::outfile->Printf("     EL  ESP   = %13.6f\n", e_el_mol );
        psi::outfile->Printf("     -------------------------------\n");
        psi::outfile->Printf("     ESP ESP   = %13.6f\n", qq       );
        psi::outfile->Printf("     -------------------------------\n");
        psi::outfile->Printf("     TOTAL VD  = %13.6f\n", e_mol_nuc + e_mol_el);
        psi::outfile->Printf("     TOTAL DV  = %13.6f\n", e_nuc_mol + e_el_mol);
        psi::outfile->Printf("     TOTAL SYM = %13.6f\n", e                   );
        psi::outfile->Printf("\n");
     }

  } 
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for electrostatic energy calculations!\n");
  }
  return e;

}
double ElectrostaticEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  double e_nuc_nuc = 0.0;
  double e_nuc_el  = 0.0; 
  double e_el_el   = 0.0;

  if (method == "DEFAULT" || method == "AO_EXPANDED") {
                            
     psi::timer_on("SOLVER: Electrostatic Energy Calculations (BENCHMARK)");
                                         
     SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
     SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
                                                                     
     // ===> Nuc(A) --- Nuc(B) <=== //
     e_nuc_nuc = wfn_union_->nuclear_repulsion_interaction_energy();
                                                                    
     // ===> Nuc(A) --- Nuc(B)  AND   Nuc(A) --- Nuc(B) <=== //
     int nbf_1 = wfn_union_->l_primary(0)->nbf();
     int nbf_2 = wfn_union_->l_primary(1)->nbf();

     std::shared_ptr<psi::Matrix> V1 = std::make_shared<psi::Matrix>("V1", nbf_1, nbf_1);
     std::shared_ptr<psi::Matrix> V2 = std::make_shared<psi::Matrix>("V2", nbf_2, nbf_2);

     psi::IntegralFactory fact_1(wfn_union_->l_primary(0), wfn_union_->l_primary(0));
     psi::IntegralFactory fact_2(wfn_union_->l_primary(1), wfn_union_->l_primary(1));

     std::shared_ptr<psi::OneBodyAOInt> oneInt;
     std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_1.spherical_transform(),
                                                                                       wfn_union_->l_primary(0),
                                                                                       wfn_union_->l_primary(0));
     std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_2.spherical_transform(),
                                                                                       wfn_union_->l_primary(1),
                                                                                       wfn_union_->l_primary(1));

     std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
     std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

     potInt_1->set_charge_field(Zxyz_2);
     potInt_2->set_charge_field(Zxyz_1);

     oneInt = potInt_1;
     oneInt->compute(V1);
     oneInt = potInt_2;
     oneInt->compute(V2);

     e_nuc_el += wfn_union_->l_wfn(0)->Da()->vector_dot(V1);
     e_nuc_el += wfn_union_->l_wfn(1)->Da()->vector_dot(V2);
     e_nuc_el += wfn_union_->l_wfn(0)->Db()->vector_dot(V1);
     e_nuc_el += wfn_union_->l_wfn(1)->Db()->vector_dot(V2);


     // ===> Nuc(A) --- Nuc(B) <=== //

     double** Da1p  = wfn_union_->l_wfn(0)->Da()->pointer();
     double** Db1p  = wfn_union_->l_wfn(0)->Db()->pointer();
     double** Da2p  = wfn_union_->l_wfn(1)->Da()->pointer();
     double** Db2p  = wfn_union_->l_wfn(1)->Db()->pointer();

     //std::shared_ptr<psi::IntegralFactory> ints_12 = std::make_shared<psi::IntegralFactory>
     //                                             (wfn_union_->l_primary(0), wfn_union_->l_primary(0),
     //                                              wfn_union_->l_primary(1), wfn_union_->l_primary(1));
     
     //std::shared_ptr<psi::TwoBodyAOInt> tei(ints_12->eri());
     //const double * buffer = tei->buffer();

     //oepdev::AllAOShellCombinationsIterator shellIter(ints_12);
     //int i, j, k, l;
     //double integral;
     //for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
     //{
     //     shellIter.compute_shell(tei);
     //     oepdev::AllAOIntegralsIterator intsIter(shellIter);
     //     for (intsIter.first(); intsIter.is_done() == false; intsIter.next())
     //     {
     //          // Grab (ij|kl) integrals and indices here
     //          i = intsIter.i();
     //          j = intsIter.j();
     //          k = intsIter.k();
     //          l = intsIter.l();
     //          integral = buffer[intsIter.index()];

     //          e_el_el += (Da1p[i][j] + Db1p[i][j]) * (Da2p[k][l] + Db2p[k][l]) * integral;
     //     }
     //}

     std::shared_ptr<psi::IntegralFactory> ints = std::make_shared<psi::IntegralFactory>
                                                  (wfn_union_->primary(), wfn_union_->primary(),
                                                   wfn_union_->primary(), wfn_union_->primary());
     
     std::shared_ptr<psi::TwoBodyAOInt> tei(ints->eri());
     const double * buffer = tei->buffer();

     oepdev::AllAOShellCombinationsIterator shellIter(ints);
     int i, j, k, l;
     double integral;

     for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
     {
          shellIter.compute_shell(tei);
          oepdev::AllAOIntegralsIterator intsIter(shellIter);
          for (intsIter.first(); intsIter.is_done() == false; intsIter.next())
          {
               i = intsIter.i();                    
               j = intsIter.j();
               k = intsIter.k();
               l = intsIter.l();

               if (i < nbf_1) {
               if (j < nbf_1) {
               if (k >= nbf_1) {
               if (l >= nbf_1) {
                   integral = buffer[intsIter.index()];
                   e_el_el += (Da1p[i][j] + Db1p[i][j]) * 
                              (Da2p[k - nbf_1][l - nbf_1] + Db2p[k - nbf_1][l - nbf_1]) * integral;
               }
               }
               }
               }

          }
     }
     // Finish //
     e = e_nuc_nuc + e_nuc_el + e_el_el;

     psi::timer_off("SOLVER: Electrostatic Energy Calculations (BENCHMARK)");

     // Print
     if (wfn_union_->options().get_int("PRINT") > 0) {
        psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
        psi::outfile->Printf("  ==>         Benchmark (DEFAULT)               <==\n\n");
        psi::outfile->Printf("     NUC NUC   = %13.6f\n", e_nuc_nuc);
        psi::outfile->Printf("     NUC EL    = %13.6f\n", e_nuc_el );
        psi::outfile->Printf("     EL  EL    = %13.6f\n", e_el_el  );
        psi::outfile->Printf("     -------------------------------\n");
        psi::outfile->Printf("     TOTAL     = %13.6f\n", e        );
        psi::outfile->Printf("\n");
     }
  } 
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for electrostatic energy calculations!\n");
  }
  return e;
}
// ===> Repulsion Energy <=== //
RepulsionEnergySolver::RepulsionEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  // Benchmarks
  methods_benchmark_.push_back("HAYES_STONE"     );
  methods_benchmark_.push_back("MURRELL_ETAL"    );
  methods_benchmark_.push_back("EFP2"            );
  // OEP-based
  methods_oepBased_ .push_back("MURRELL_ETAL_MIX");
  methods_oepBased_ .push_back("MURRELL_ETAL_ESP");
}
RepulsionEnergySolver::~RepulsionEnergySolver() 
{

}
double RepulsionEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  if      (method == "DEFAULT" || 
           method == "MURRELL_ETAL_MIX")  e = compute_oep_based_murrell_etal_mix();
  else if (method == "MURRELL_ETAL_ESP")  e = compute_oep_based_murrell_etal_esp();
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
           method == "HAYES_STONE" ) e = compute_benchmark_hayes_stone();
  else if (method == "MURRELL_ETAL") e = compute_benchmark_murrell_etal();
  else if (method == "EFP2")         e = compute_benchmark_efp2();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for repulsion energy calculations!\n");
  }
  return e;
}
double RepulsionEnergySolver::compute_benchmark_hayes_stone() {

  // ===> Start computations <=== //
  psi::timer_on ("SOLVER: Repulsion Energy Calculations (Hayes-Stone (1984))");

  double e_1, e_2, e_ex;

  // ===> One electron part <=== //

  int nbf   = wfn_union_->basisset()->nbf();

  // ---> Allocate <--- //
  std::shared_ptr<psi::Matrix> VaoA      = std::make_shared<psi::Matrix>("VaoA" , nbf, nbf);
  std::shared_ptr<psi::Matrix> VaoB      = std::make_shared<psi::Matrix>("VaoB" , nbf, nbf);
  std::shared_ptr<psi::Matrix> Sao       = std::make_shared<psi::Matrix>("Sao"  , nbf, nbf);
  std::shared_ptr<psi::Matrix> Tao       = std::make_shared<psi::Matrix>("Tao"  , nbf, nbf);
  psi::IntegralFactory fact(wfn_union_->basisset(), wfn_union_->basisset());

  //std::shared_ptr<psi::Matrix> VaoA21    = std::make_shared<psi::Matrix>("VaoA(2,1)" , nbf_2, nbf_1);
  //std::shared_ptr<psi::Matrix> VaoB12    = std::make_shared<psi::Matrix>("VaoB(1,2)" , nbf_1, nbf_2);
  //std::shared_ptr<psi::Matrix> Sao12     = std::make_shared<psi::Matrix>("Sao(1,2)"  , nbf_1, nbf_2);
  //std::shared_ptr<psi::Matrix> Tao12     = std::make_shared<psi::Matrix>("Tao(1,2)"  , nbf_1, nbf_2);

  psi::IntegralFactory fact_12(wfn_union_->l_primary(0), wfn_union_->l_primary(1));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1), wfn_union_->l_primary(0));

  std::shared_ptr<psi::OneBodyAOInt> oneInt;
  std::shared_ptr<psi::PotentialInt> potInt_1 = std::make_shared<psi::PotentialInt>(fact_12.spherical_transform(),
                                                                                    wfn_union_->l_primary(0),
                                                                                    wfn_union_->l_primary(1));
  std::shared_ptr<psi::PotentialInt> potInt_2 = std::make_shared<psi::PotentialInt>(fact_21.spherical_transform(),
                                                                                    wfn_union_->l_primary(1),
                                                                                    wfn_union_->l_primary(0));
  std::shared_ptr<psi::PotentialInt> potInt_1n= std::make_shared<psi::PotentialInt>(fact.spherical_transform(),
                                                                                    wfn_union_->basisset(),
                                                                                    wfn_union_->basisset());
  std::shared_ptr<psi::PotentialInt> potInt_2n= std::make_shared<psi::PotentialInt>(fact.spherical_transform(),
                                                                                    wfn_union_->basisset(),
                                                                                    wfn_union_->basisset());

  std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact.ao_overlap());
  std::shared_ptr<psi::OneBodyAOInt> kinInt(fact.ao_kinetic());

  std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
  std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

  potInt_1n->set_charge_field(Zxyz_2);
  potInt_2n->set_charge_field(Zxyz_1);

  // ---> Compute one-electron integrals <--- //
  oneInt = potInt_1n;
  oneInt->compute(VaoB);
  oneInt = potInt_2n;
  oneInt->compute(VaoA);
  ovlInt->compute(Sao);
  kinInt->compute(Tao);

  // ---> Accumulate one electron part <--- //
  Tao->add(VaoA);
  Tao->add(VaoB);
  Tao->scale(2.0);

  // ---> Transform one electron and overlap matrices to MO basis <--- //
  std::shared_ptr<psi::Matrix> Smo = psi::Matrix::triplet(wfn_union_->Ca_subset("AO","OCC"), 
                                                          Sao, 
                                                          wfn_union_->Ca_subset("AO","OCC"),
                                                          true, false, false);
  std::shared_ptr<psi::Matrix> Tmo = psi::Matrix::triplet(wfn_union_->Ca_subset("AO","OCC"), 
                                                          Tao, 
                                                          wfn_union_->Ca_subset("AO","OCC"),
                                                          true, false, false);

  // ---> Invert the overlap matrix in MO basis <--- //
  Smo->invert();

  // ---> Finalize with one-electron term <--- //  
  std::shared_ptr<psi::Matrix> Imo = std::make_shared<psi::Matrix>(Smo);
  Imo->identity();
  Smo->subtract(Imo);
  e_1 = Tmo->vector_dot(Smo);

  // ===> Two electron part <=== //
  Smo->add(Imo);
  e_2 = 0.0;
  // ---> Loop over ERI's in MO space <--- //

  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  dpdbuf4 buf;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  psio->tocprint(PSIF_LIBTRANS_DPD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[O,O]"  ), integrals->DPD_ID("[O,O]"  ),
                          integrals->DPD_ID("[O>=O]+"), integrals->DPD_ID("[O>=O]+"), 0, "MO Ints (OO|OO)");

  double integral, v;
  double** t = Smo->pointer();
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
       for (int kl = 0; kl < buf.params->rowtot[h]; ++kl) {
            int k = buf.params->roworb[h][kl][0];
            int l = buf.params->roworb[h][kl][1];
            for (int mn = 0; mn < buf.params->coltot[h]; ++mn) {
                 int m = buf.params->colorb[h][mn][0];
                 int n = buf.params->colorb[h][mn][1];
                 //psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);

                 integral = buf.matrix[h][kl][mn];

                 v = 2.0*t[k][l]*t[m][n] - t[k][n]*t[m][l];
                 if (k != l) v+= 2.0*t[k][l]*t[m][n] - t[l][n]*t[k][m];
                 if (m != n) v+= 2.0*t[k][l]*t[m][n] - t[l][m]*t[k][n];
                 if ((k != l) && (m != n)) v+= 2.0*t[k][l]*t[m][n] - t[k][m]*t[l][n];

                 if ((k == l) && (m == n)) v -= 2.0;
                 if ((k == n) && (l == m)) v += 2.0;

                 e_2 += integral * v;
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);

  psi::timer_off("SOLVER: Repulsion Energy Calculations (Hayes-Stone (1984))");

  // ===> Compute Exchange Energy <=== //
  psi::timer_on ("SOLVER: HF Exchange Energy Calculations");

  dpdbuf4 buf_1212;

  global_dpd_->buf4_init(&buf_1212, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[1,2]"  ),
                          integrals->DPD_ID("[1,2]"  ), integrals->DPD_ID("[1,2]"  ), 0, "MO Ints (12|12)");

  e_ex = 0.0;
  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf_1212, h);
       global_dpd_->buf4_mat_irrep_rd(&buf_1212, h);
       for (int pq = 0; pq < buf_1212.params->rowtot[h]; ++pq) {
            int p = buf_1212.params->roworb[h][pq][0];
            int q = buf_1212.params->roworb[h][pq][1];
            for (int rs = 0; rs < buf_1212.params->coltot[h]; ++rs) {
                 int r = buf_1212.params->colorb[h][rs][0];
                 int s = buf_1212.params->colorb[h][rs][1];
                 e_ex += buf_1212.matrix[h][pq][rs];
                 //psi::outfile->Printf("(%2d %2d | %2d %2d) = %16.10f\n", p, q, r, s, buf_1212.matrix[h][pq][rs]);
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf_1212, h);
  }
  global_dpd_->buf4_close(&buf_1212);
  e_ex *= -2.0;

  psi::timer_off("SOLVER: HF Exchange Energy Calculations");

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Exchange-Repulsion energy calculations <==\n");
     psi::outfile->Printf("  ==>         Benchmark (Hayes-Stone)                <==\n\n");
     psi::outfile->Printf("     E_REP_1   = %13.6f\n", e_1      );
     psi::outfile->Printf("     E_REP_2   = %13.6f\n", e_2      );
     psi::outfile->Printf("     E_REP     = %13.6f\n", e_1+e_2  );
     psi::outfile->Printf("     E_EX      = %13.6f\n", e_ex     );
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("     E_EXREP   = %13.6f\n", e_1+e_2+e_ex);
     psi::outfile->Printf("\n");
  }

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  return e_1+e_2+e_ex;
}
double RepulsionEnergySolver::compute_benchmark_murrell_etal() {
  psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell et al.)");
  // ===> One electron part <=== //
  // ===> One electron part <=== //
  psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell et al.)");
  throw psi::PSIEXCEPTION("ERROR: MURRELL_ETAL is not yet implemented!\n");
}
double RepulsionEnergySolver::compute_benchmark_efp2() {
  psi::timer_on ("SOLVER: Repulsion Energy Calculations (EFP2)");
  psi::timer_off("SOLVER: Repulsion Energy Calculations (EFP2)");
  throw psi::PSIEXCEPTION("ERROR: EFP2 is not yet implemented!\n");
}
double RepulsionEnergySolver::compute_oep_based_murrell_etal_mix() {
  psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell-OEP:S1-DF/S2-ESP)");
  psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell-OEP:S1-DF/S2-ESP)");
  throw psi::PSIEXCEPTION("ERROR: MURRELL_ETAL_MIX is not yet implemented!\n");
}
double RepulsionEnergySolver::compute_oep_based_murrell_etal_esp() {
  psi::timer_on ("SOLVER: Repulsion Energy Calculations (Murrell-OEP:ESP)");
  psi::timer_off("SOLVER: Repulsion Energy Calculations (Murrell-OEP:ESP)");
  throw psi::PSIEXCEPTION("ERROR: MURRELL_ETAL_ESP is not yet implemented!\n");
}
// Build: factory static method 
std::shared_ptr<OEPDevSolver> OEPDevSolver::build(const std::string& target, SharedWavefunctionUnion wfn_union)
{
   std::shared_ptr<OEPDevSolver> solver;
   if      (target == "ELECTROSTATIC ENERGY") solver = std::make_shared<ElectrostaticEnergySolver>(wfn_union);
   else if (target == "REPULSION ENERGY"    ) solver = std::make_shared<    RepulsionEnergySolver>(wfn_union);
   else throw psi::PSIEXCEPTION("OEPDEV: Error. OEPDevSolver: Incorrect targer property chosen!\n");
   return solver;
}

