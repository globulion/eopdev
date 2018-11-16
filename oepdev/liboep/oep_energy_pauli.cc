#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

// <============== Repulsion Energy ==============> //

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, 
                                                       SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::~RepulsionEnergyOEPotential() 
{
  delete[] vec_otto_ladik_s2_;
}
void RepulsionEnergyOEPotential::common_init() 
{
   int n1 = wfn_->Ca_subset("AO","OCC")->ncol();
   int n2 = auxiliary_->nbf();
   int n3 = wfn_->molecule()->natom();

   SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G(S^{-1})", n2, n1);
   SharedMatrix mat_2 = std::make_shared<psi::Matrix>("G(S^{-2})", n3, n1);

   OEPType type_1 = {"Murrell-etal.S1", true , n1, mat_1};
   OEPType type_2 = {"Otto-Ladik.S2"  , false, n1, mat_2};

   oepTypes_[type_1.name] = type_1; 
   oepTypes_[type_2.name] = type_2;

   //
   vec_otto_ladik_s2_ = new double[n1];
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.S1") compute_murrell_etal_s1();
  else if (oepType ==   "Otto-Ladik.S2") compute_otto_ladik_s2();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void RepulsionEnergyOEPotential::compute_murrell_etal_s1() 
{
   //psi::timer_on("OEP    E(Paul) Murrell-etal S1  ");

   // ---> Determine the target basis set for generalized density fitting <--- //
   std::shared_ptr<psi::BasisSet> target = intermediate_;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") target = auxiliary_;

   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> Ca_occ = wfn_->Ca_subset("AO","OCC");

   psi::IntegralFactory fact_1(primary_, target, primary_, target);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, target);
  
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());

   // ===> Compute One-Electron Integrals <=== //
   potInt->compute(Vao);

   // ===> Compute Nuclear Contribution to DF Vector <=== //
   std::shared_ptr<psi::Matrix> V = psi::Matrix::doublet(Vao, Ca_occ, true, false);

   // ===> Compute Electronic Contribution to DF Vector <=== //
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact_2.eri());                                           
   const double * buffer = tei->buffer();
                                                                                                  
   std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact_2, "ALL");
   int i, j, k, l;
   double integral, dij, dik;
   
   double** v = V->pointer();  
   double** c = Ca_occ->pointer();
   double** d = wfn_->Da()->pointer();
   for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
   {
        shellIter->compute_shell(tei);
        std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
        for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
        {
             i = intsIter->i();  // \mu      : n
             j = intsIter->j();  // \nu      : n
             k = intsIter->k();  // \alpha   : n
             l = intsIter->l();  // \xi      : Q

             integral = buffer[intsIter->index()];

             dij = d[i][j]; dik = d[i][k];
             for (int a = 0; a < wfn_->doccpi()[0]; ++a) 
             {
                v[l][a] += (2.0 * c[k][a]*dij - c[j][a]*dik) * integral;
             }
        }
   }

   // ===> Perform The Generalized Density Fitting <=== // 
   std::shared_ptr<oepdev::GeneralizedDensityFit> gdf;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, V);}
   else                                                {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, intermediate_, V);}
   std::shared_ptr<psi::Matrix> G = gdf->compute();
   
   // ===> Save and Finish <=== //
   oepTypes_.at("Murrell-etal.S1").matrix->copy(G);
   if (options_.get_int("PRINT") > 1) G->print();
   //psi::timer_off("OEP    E(Paul) Murrell-etal S1  ");
}
void RepulsionEnergyOEPotential::compute_otto_ladik_s2() 
{
      //psi::timer_on("OEP    E(Paul) Otto-Ladik S2    ");
      std::shared_ptr<OEPotential3D<OEPotential>> oeps3d = this->make_oeps3d("Otto-Ladik.S2");
      oeps3d->compute();
      ESPSolver esp(oeps3d);
      esp.set_charge_sums(0.5);
      esp.compute();
     
      for (int i=0; i<esp.charges()->nrow(); ++i) {
      for (int o=0; o<oepTypes_["Otto-Ladik.S2"].n; ++o) {
           oepTypes_["Otto-Ladik.S2"].matrix->set(i, o, esp.charges()->get(i, o));
      }}
      if (options_.get_int("PRINT") > 1) esp.charges()->print();
      //psi::timer_off("OEP    E(Paul) Otto-Ladik S2    ");
}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) 
{
   if (oepType == "Otto-Ladik.S2") {
       this->compute_3D_otto_ladik_s2(x, y, z);
       // Assign final value
       for (int o = 0; o < oepTypes_["Otto-Ladik.S2"].n; ++o) v->set(o, vec_otto_ladik_s2_[o]);
   }
   else if (oepType == "Murrell-etal.S1" ) {/* nothing to do here */}
   else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
   }

}
void RepulsionEnergyOEPotential::compute_3D_otto_ladik_s2(const double& x, const double& y, const double& z)
{
  double val = 0.0;
  // ===> Nuclear contribution <=== //
  for (int i=0; i<wfn_->molecule()->natom(); ++i) {
       val+= (double)wfn_->molecule()->Z(i) /
                 sqrt( pow(wfn_->molecule()->x(i)-x, 2.0) 
                     + pow(wfn_->molecule()->y(i)-y, 2.0) 
                     + pow(wfn_->molecule()->z(i)-z, 2.0) );
  }

  // ===> Electronic contribution <=== //
  potInt_->set_charge_field(x, y, z);
  OEInt_ = potInt_;
  OEInt_->compute(potMat_);
  std::shared_ptr<psi::Matrix> potMO = psi::Matrix::triplet(cOcc_, potMat_, cOcc_, true, false, false);
  potMat_->zero();
  //for (int o=0; o<oepTypes_["Otto-Ladik.S2"].n; ++o) val += 2.0 * potMO->get(o, o);
  val += 2.0 * potMO->trace();

  for (int o=0; o<oepTypes_["Otto-Ladik.S2"].n; ++o) {
     vec_otto_ladik_s2_[o] = val - 0.5 * potMO->get(o, o);
  }
}
void RepulsionEnergyOEPotential::print_header(void) const {}
