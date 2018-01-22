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
  methods_oepBased_.push_back("DEFAULT");
  methods_benchmark_.push_back("DEFAULT");
}

ElectrostaticEnergySolver::~ElectrostaticEnergySolver() 
{

}

double ElectrostaticEnergySolver::compute_oep_based(const std::string& method) 
{
  double e = 0.0;
  double e_nuc_nuc = 0.0;
  double e_mol_el  = 0.0;

  if (method == "DEFAULT") {

     psi::timer_on("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");

     SharedWavefunction wfn_1 = wfn_union_->l_wfn(0);
     SharedWavefunction wfn_2 = wfn_union_->l_wfn(1);
     SharedOEPotential oep_1 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_1, wfn_union_->options());
     SharedOEPotential oep_2 = oepdev::OEPotential::build("ELECTROSTATIC ENERGY", wfn_2, wfn_union_->options());
                                                                     
     // ===> Nuc(A) --- Nuc(B) <=== //
     e_nuc_nuc = wfn_union_->nuclear_repulsion_interaction_energy();

     // ===> [Nuc+El](A) --- El(B) <=== //
     
     // Finish
     e = e_nuc_nuc + e_mol_el;

     psi::timer_off("SOLVER: Electrostatic Energy Calculations (OEP-BASED)");

     // Print
     if (wfn_union_->options().get_int("PRINT") > 0) {
        psi::outfile->Printf("  ==> SOLVER: Enectrostatic energy calculations <==\n");
        psi::outfile->Printf("  ==>         OEP-based Model (DEFAULT)         <==\n\n");
        psi::outfile->Printf("     NUC NUC   = %13.6f\n", e_nuc_nuc);
        psi::outfile->Printf("     MOL EL    = %13.6f\n", e_mol_el );
        psi::outfile->Printf("     -------------------------------\n");
        psi::outfile->Printf("     TOTAL     = %13.6f\n", e        );
        psi::outfile->Printf("\n");
     }

  } 
  else 
  {
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for electrostatic calculations!\n");
  }
  return e;

}
double ElectrostaticEnergySolver::compute_benchmark(const std::string& method) 
{
  double e = 0.0;
  double e_nuc_nuc = 0.0;
  double e_nuc_el  = 0.0; 
  double e_el_el   = 0.0;

  if (method == "DEFAULT") {
                            
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
     throw psi::PSIEXCEPTION("Error. Incorrect benchmark method specified for electrostatic calculations!\n");
  }
  return e;
}



std::shared_ptr<OEPDevSolver> OEPDevSolver::build(const std::string& target, SharedWavefunctionUnion wfn_union)
{
   std::shared_ptr<OEPDevSolver> solver;
   if     (target == "ELECTROSTATIC ENERGY") solver = std::make_shared<ElectrostaticEnergySolver>(wfn_union);
   else throw psi::PSIEXCEPTION("OEPDEV: Error. OEPDevSolver: Incorrect targer property chosen!\n");
   return solver;
}

