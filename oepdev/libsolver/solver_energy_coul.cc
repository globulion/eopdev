#include "solver.h"
#include "../lib3d/dmtp.h"
#include "psi4/libpsi4util/process.h"

using namespace std;
using namespace psi;
using namespace oepdev;

ElectrostaticEnergySolver::ElectrostaticEnergySolver(SharedWavefunctionUnion wfn_union)
 : OEPDevSolver(wfn_union)
{
  methods_oepBased_ .push_back("ESP_SYMMETRIZED");
  methods_benchmark_.push_back("AO_EXPANDED"    );
  methods_benchmark_.push_back("MO_EXPANDED"    );
}
ElectrostaticEnergySolver::~ElectrostaticEnergySolver() {}
double ElectrostaticEnergySolver::compute_oep_based_camm(){

  // ===> [Nuc+El](A) --- El(B) <=== //
  SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
  SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
  SharedOEPotential oep_1 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_1, wfn_union_->options());
  SharedOEPotential oep_2 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_2, wfn_union_->options());

  oep_1->compute();
  oep_2->compute();

  psi::timer_on("Solver E(Coul) OEP-Based:CAMM   ");

  double e1 = (oep_1->oep("V").dmtp->energy(oep_2->oep("V").dmtp))->level(oepdev::MultipoleConvergence::R1)->get(0,0);
  double e2 = (oep_1->oep("V").dmtp->energy(oep_2->oep("V").dmtp))->level(oepdev::MultipoleConvergence::R2)->get(0,0);
  double e3 = (oep_1->oep("V").dmtp->energy(oep_2->oep("V").dmtp))->level(oepdev::MultipoleConvergence::R3)->get(0,0);
  double e4 = (oep_1->oep("V").dmtp->energy(oep_2->oep("V").dmtp))->level(oepdev::MultipoleConvergence::R4)->get(0,0);
  double e5 = (oep_1->oep("V").dmtp->energy(oep_2->oep("V").dmtp))->level(oepdev::MultipoleConvergence::R5)->get(0,0);

  double energy = e5;

  psi::timer_off("Solver E(Coul) OEP-Based:CAMM   ");

  // Save
  psi::Process::environment.globals["EINT COUL CAMM R-1"] = e1;
  psi::Process::environment.globals["EINT COUL CAMM R-2"] = e2;
  psi::Process::environment.globals["EINT COUL CAMM R-3"] = e3;
  psi::Process::environment.globals["EINT COUL CAMM R-4"] = e4;
  psi::Process::environment.globals["EINT COUL CAMM R-5"] = e5;


  // Print
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
     psi::outfile->Printf("  ==>         OEP-based Model (CAMM)            <==\n\n");
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("     E(R-1)    = %13.6f\n", e1  );
     psi::outfile->Printf("     E(R-2)    = %13.6f\n", e2  );
     psi::outfile->Printf("     E(R-3)    = %13.6f\n", e3  );
     psi::outfile->Printf("     E(R-4)    = %13.6f\n", e4  );
     psi::outfile->Printf("     E(R-5)    = %13.6f\n", e5  );
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("     E         = %13.6f\n", energy   );
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("\n");
  }
return energy;
}

double ElectrostaticEnergySolver::compute_oep_based_esp_symmetrized(){

  double e = 0.0;
  double e_nuc_mol = 0.0;
  double e_mol_nuc = 0.0;
  double e_el_mol  = 0.0;
  double e_mol_el  = 0.0;
  double qq        = 0.0;

  // ===> [Nuc+El](A) --- El(B) <=== //
  SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
  SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
  SharedOEPotential oep_1 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_1, wfn_union_->options());
  SharedOEPotential oep_2 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_2, wfn_union_->options());

  oep_1->compute();
  oep_2->compute();

  oep_1->print_header();
  oep_2->print_header();

  //psi::timer_on("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");
  psi::timer_on("Solver E(Coul) OEP-Based:MIX    ");

  int nbf_1 = wfn_union_->l_primary(0)->nbf();
  int nbf_2 = wfn_union_->l_primary(1)->nbf();

  std::shared_ptr<psi::Matrix> V1 = std::make_shared<psi::Matrix>("V1", nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> V2 = std::make_shared<psi::Matrix>("V2", nbf_2, nbf_2);

  psi::IntegralFactory fact_1(wfn_union_->l_primary(0));
  psi::IntegralFactory fact_2(wfn_union_->l_primary(1));

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
       e_mol_nuc += Qxyz_1->get(i, 0) * (double)oep_2->wfn()->molecule()->Z(j) /
            sqrt( pow( Qxyz_1->get(i,1) - oep_2->wfn()->molecule()->x(j), 2.0) +
                  pow( Qxyz_1->get(i,2) - oep_2->wfn()->molecule()->y(j), 2.0) +
                  pow( Qxyz_1->get(i,3) - oep_2->wfn()->molecule()->z(j), 2.0) );
  }
  }
  // ===> Nuc(A) --- ESP(B) <=== //
  for (int i=0; i < oep_1->wfn()->molecule()->natom(); ++i) {
  for (int j=0; j < Qxyz_2->nrow(); ++j) {
       e_nuc_mol += Qxyz_2->get(j, 0) * (double)oep_1->wfn()->molecule()->Z(i) /
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

  //psi::timer_off("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");
  psi::timer_off("Solver E(Coul) OEP-Based:MIX    ");

  // ===> ESP(A) ---- ESP(B) <=== //
  //psi::timer_on ("SOLVER: Electrostatic Energy Calculations (ESP-ESP)");
  psi::timer_on("Solver E(Coul) ESP              ");

  for (int i=0; i < Qxyz_1->nrow(); ++i) {
  for (int j=0; j < Qxyz_2->nrow(); ++j) {
       qq += Qxyz_1->get(i, 0) * Qxyz_2->get(j, 0) /
            sqrt( pow( Qxyz_1->get(i,1) - Qxyz_2->get(j,1), 2.0) +
                  pow( Qxyz_1->get(i,2) - Qxyz_2->get(j,2), 2.0) +
                  pow( Qxyz_1->get(i,3) - Qxyz_2->get(j,3), 2.0) );
  }
  }
  //psi::timer_off("SOLVER: Electrostatic Energy Calculations (ESP-ESP)");
  psi::timer_off("Solver E(Coul) ESP              ");

  // Save
  psi::Process::environment.globals["EINT COUL ESP"] = e;

  // Print
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
     psi::outfile->Printf("  ==>         OEP-based Model (ESP-Symmetrized) <==\n\n");
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
return e;
}
double ElectrostaticEnergySolver::compute_oep_based(const std::string& method) 
{
  double e;
  if      (method == "DEFAULT" 
	|| method == "ESP_SYMMETRIZED") {e = compute_oep_based_esp_symmetrized();}
  else if (method == "CAMM")            {e = compute_oep_based_camm();}
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect OEP-based method specified for electrostatic energy calculations!\n");
  }
  return e;

}
double ElectrostaticEnergySolver::compute_benchmark_mo_expanded(){
  double e = 0.0;
  double e_nuc_nuc = wfn_union_->nuclear_repulsion_interaction_energy();
  double e_nuc_el  = 0.0; 
  double e_el_el   = 0.0;

  //psi::timer_on ("SOLVER: Electrostatic Energy Calculations (MO-expanded)");
  psi::timer_on("Solver E(Coul) MO-Expanded      ");

  // ===> One electron part <=== //

  int nbf   = wfn_union_->basisset()->nbf();

  // ---> Allocate <--- //
  std::shared_ptr<psi::Matrix> VaoA      = std::make_shared<psi::Matrix>("VaoA" , nbf, nbf);
  std::shared_ptr<psi::Matrix> VaoB      = std::make_shared<psi::Matrix>("VaoB" , nbf, nbf);
  psi::IntegralFactory fact(wfn_union_->basisset());

  psi::IntegralFactory fact_12(wfn_union_->l_primary(0));
  psi::IntegralFactory fact_21(wfn_union_->l_primary(1));

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

  std::shared_ptr<psi::Matrix> Zxyz_1 = std::make_shared<psi::Matrix>(potInt_1->charge_field());
  std::shared_ptr<psi::Matrix> Zxyz_2 = std::make_shared<psi::Matrix>(potInt_2->charge_field());

  potInt_1n->set_charge_field(Zxyz_2);
  potInt_2n->set_charge_field(Zxyz_1);

  // ---> Compute one-electron integrals <--- //
  oneInt = potInt_1n;
  oneInt->compute(VaoB);
  oneInt = potInt_2n;
  oneInt->compute(VaoA);

  // ---> Transform one electron contributions to MO basis <--- //
  std::shared_ptr<psi::Matrix> Ca_occ = wfn_union_->Ca_subset("AO","OCC");
  std::shared_ptr<psi::Matrix> Ca_A   = std::make_shared<psi::Matrix>("Ca_A", nbf, wfn_union_->l_ndocc(0));
  std::shared_ptr<psi::Matrix> Ca_B   = std::make_shared<psi::Matrix>("Ca_B", nbf, wfn_union_->l_ndocc(1));
  for (int i=0; i<nbf; ++i) {
       for (int ja=0; ja<wfn_union_->l_ndocc(0); ++ja) {
            Ca_A->set(i, ja, Ca_occ->get(i, ja));
       }
       for (int jb=0; jb<wfn_union_->l_ndocc(1); ++jb) {
            Ca_B->set(i, jb, Ca_occ->get(i, jb+wfn_union_->l_ndocc(0)));
       }
  }
  std::shared_ptr<psi::Matrix> VmoA = psi::Matrix::triplet(Ca_B, VaoA, Ca_B, true, false, false);
  std::shared_ptr<psi::Matrix> VmoB = psi::Matrix::triplet(Ca_A, VaoB, Ca_A, true, false, false);

  // ---> Finalize with one-electron term <--- //  
  e_nuc_el = 2.0 * (VmoA->trace() + VmoB->trace());


  // ===> Two electron part <=== //

  // ---> Loop over ERI's in MO space <--- //
  std::shared_ptr<psi::IntegralTransform> integrals = wfn_union_->integrals(); 
  dpd_set_default(integrals->get_dpd_id());
  dpdbuf4 buf;
  std::shared_ptr<psi::PSIO> psio = psi::PSIO::shared_object();
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  global_dpd_->buf4_init(&buf, PSIF_LIBTRANS_DPD, 0, 
                          integrals->DPD_ID("[I,I]"  ), integrals->DPD_ID("[J,J]"  ),
                          integrals->DPD_ID("[I>=I]+"), integrals->DPD_ID("[J>=J]+"), 0, "MO Ints (II|JJ)");

  for (int h = 0; h < wfn_union_->nirrep(); ++h) {
       global_dpd_->buf4_mat_irrep_init(&buf, h);
       global_dpd_->buf4_mat_irrep_rd(&buf, h);
       for (int kl = 0; kl < buf.params->rowtot[h]; ++kl) {
            int k = buf.params->roworb[h][kl][0];
            int l = buf.params->roworb[h][kl][1];
            for (int mn = 0; mn < buf.params->coltot[h]; ++mn) {
                 int m = buf.params->colorb[h][mn][0];
                 int n = buf.params->colorb[h][mn][1];
                 //psi::outfile->Printf(" Yint: (%2d %2d | %2d %2d) = %16.10f\n", k, l, m, n, buf.matrix[h][kl][mn]);
                 if ((k==l) && (m==n)) e_el_el += buf.matrix[h][kl][mn];
            }
       }
       global_dpd_->buf4_mat_irrep_close(&buf, h);
  }
  global_dpd_->buf4_close(&buf);
  e_el_el *= 4.0;

  // ---> Close the DPD file <--- //
  psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

  //psi::timer_off("SOLVER: Electrostatic Energy Calculations (MO-expanded)");
  psi::timer_off("Solver E(Coul) MO-Expanded      ");

  // ---> Sum <--- //
  e = e_nuc_nuc + e_nuc_el + e_el_el;

  // ---> Save <--- //
  psi::Process::environment.globals["EINT COUL EXACT"] = e;

  // ---> Print <--- //
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Electrostatic energy calculations <==\n");
     psi::outfile->Printf("  ==>         Benchmark (MO-Expanded)           <==\n\n");
     psi::outfile->Printf("     NUC NUC   = %13.6f\n", e_nuc_nuc);
     psi::outfile->Printf("     NUC EL    = %13.6f\n", e_nuc_el );
     psi::outfile->Printf("     EL  EL    = %13.6f\n", e_el_el  );
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("     E_TOTAL   = %13.6f\n", e        );
     psi::outfile->Printf("\n");
  }


  return e;
}
double ElectrostaticEnergySolver::compute_benchmark_ao_expanded(){

  double e = 0.0;
  double e_nuc_nuc = 0.0;
  double e_nuc_el  = 0.0; 
  double e_el_el   = 0.0;

  //psi::timer_on("SOLVER: Electrostatic Energy Calculations (BENCHMARK)");
  psi::timer_on("Solver E(Coul) AO-Expanded      ");
                                      
  SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
  SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
                                                                  
  // ===> Nuc(A) --- Nuc(B) <=== //
  e_nuc_nuc = wfn_union_->nuclear_repulsion_interaction_energy();
                                                                 
  // ===> Nuc(A) --- Nuc(B)  AND   Nuc(A) --- Nuc(B) <=== //
  int nbf_1 = wfn_union_->l_primary(0)->nbf();
  int nbf_2 = wfn_union_->l_primary(1)->nbf();

  std::shared_ptr<psi::Matrix> V1 = std::make_shared<psi::Matrix>("V1", nbf_1, nbf_1);
  std::shared_ptr<psi::Matrix> V2 = std::make_shared<psi::Matrix>("V2", nbf_2, nbf_2);

  psi::IntegralFactory fact_1(wfn_union_->l_primary(0));
  psi::IntegralFactory fact_2(wfn_union_->l_primary(1));

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

  std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(ints, "ALL");
  int i, j, k, l;
  double integral;

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

            if (i < nbf_1) {
            if (j < nbf_1) {
            if (k >= nbf_1) {
            if (l >= nbf_1) {
                integral = buffer[intsIter->index()];
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

  //psi::timer_off("SOLVER: Electrostatic Energy Calculations (BENCHMARK)");
  psi::timer_off("Solver E(Coul) AO-Expanded      ");

  // ---> Save <--- //
  psi::Process::environment.globals["EINT COUL EXACT"] = e;

  // Print
  if (wfn_union_->options().get_int("PRINT") > 0) {
     psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
     psi::outfile->Printf("  ==>         Benchmark (AO-Expanded)           <==\n\n");
     psi::outfile->Printf("     NUC NUC   = %13.6f\n", e_nuc_nuc);
     psi::outfile->Printf("     NUC EL    = %13.6f\n", e_nuc_el );
     psi::outfile->Printf("     EL  EL    = %13.6f\n", e_el_el  );
     psi::outfile->Printf("     -------------------------------\n");
     psi::outfile->Printf("     TOTAL     = %13.6f\n", e        );
     psi::outfile->Printf("\n");
  }
  return e;
}
double ElectrostaticEnergySolver::compute_benchmark(const std::string& method) 
{
  double e;
  if      (method == "DEFAULT" || 
           method == "AO_EXPANDED") e = compute_benchmark_ao_expanded();
  else if (method == "MO_EXPANDED") e = compute_benchmark_mo_expanded();
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for electrostatic energy calculations!\n");
  }
  return e;
}
