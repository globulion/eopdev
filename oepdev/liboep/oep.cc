#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

OEPotential::OEPotential(SharedWavefunction wfn, Options& options) 
   : wfn_(wfn),
     options_(options)
{
    common_init();
}
OEPotential::OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options) 
   : wfn_(wfn),
     auxiliary_(auxiliary),
     intermediate_(intermediate),
     options_(options)
{
    common_init();
}
void OEPotential::common_init(void) 
{
   name_        = "default";
   primary_     = wfn_->basisset();
   intsFactory_ = std::make_shared<psi::IntegralFactory>(primary_, primary_);
   potMat_      = std::make_shared<psi::Matrix>("Potential Integrals", primary_->nbf(), primary_->nbf());
   potInt_      = std::make_shared<oepdev::PotentialInt>(intsFactory_->spherical_transform(), primary_, primary_, 0);
}
std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY"  )  oep = std::make_shared< ElectrostaticEnergyOEPotential>(wfn, options);
   else if (category == "REPULSION ENERGY"      )  oep = std::make_shared<     RepulsionEnergyOEPotential>(wfn, options);
   else if (category == "CHARGE TRANSFER ENERGY")  oep = std::make_shared<ChargeTransferEnergyOEPotential>(wfn, options);
   else if (category == "EET COUPLING CONSTANT" )  oep = std::make_shared<         EETCouplingOEPotential>(wfn, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}
std::shared_ptr<OEPotential> OEPotential::build(const std::string& category, SharedWavefunction wfn, 
                                                                           SharedBasisSet auxiliary, SharedBasisSet intermediate, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY"  )  throw PSIEXCEPTION("OEPDEV: OEPotential build. DF-based OEP not available for Electrostatic Energy!");
   else if (category == "REPULSION ENERGY"      )  oep = std::make_shared< RepulsionEnergyOEPotential>(wfn, auxiliary, intermediate, options);
   else if (category == "CHARGE TRANSFER ENERGY")  oep = std::make_shared< ChargeTransferEnergyOEPotential>(wfn, auxiliary, intermediate, options);
   else if (category == "EET COUPLING CONSTANT" )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, auxiliary, intermediate, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}
OEPotential::~OEPotential() {}
void OEPotential::compute(const std::string& oepType) {}
void OEPotential::compute(void) { for ( auto const& oepType : oepTypes_ ) this->compute(oepType.second.name); }
void OEPotential::write_cube(const std::string& oepType, const std::string& fileName) 
{
   OEPotential3D<OEPotential> oeps3d(oepTypes_[oepType].n, 60, 60, 60, 10.0, 10.0, 10.0, shared_from_this(), oepType, options_);
   oeps3d.write_cube_file(fileName);
}
void OEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void OEPotential::rotate(const Matrix& rotmat) {}
void OEPotential::translate(const Vector& trans) {}
void OEPotential::superimpose(const Matrix& refGeometry,
                              const std::vector<int>& supList,
                              const std::vector<int>& reordList) {}
void OEPotential::print_header(void) const {}


// <============== Electrostatic Energy (demo class) ==============> //

ElectrostaticEnergyOEPotential::ElectrostaticEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

ElectrostaticEnergyOEPotential::~ElectrostaticEnergyOEPotential() {}
void ElectrostaticEnergyOEPotential::common_init() 
{
   // In this demo OEPotential, there is only one OEP type. In addition,
   // only one OEP exists within this type: the electrostatic potential of entire molecule
   SharedMatrix mat = std::make_shared<psi::Matrix>("OEP V", wfn_->molecule()->natom(), 1);
   OEPType type = {"V", false, 1, mat};
   oepTypes_["V"] = type;
}

void ElectrostaticEnergyOEPotential::compute(const std::string& oepType) 
{
  if (oepType == "V" || oepType == "TOTAL") {

      psi::timer_on("OEPDEV: Electrostatic Energy OEP -> fitting ESP charges");
      SharedField3D potential = oepdev::Field3D::build("ELECTROSTATIC POTENTIAL", 
                                                   options_.get_int   ("ESP_NPOINTS_PER_ATOM") * wfn_->molecule()->natom(), 
                                                   options_.get_double("ESP_PAD_SPHERE"      ), 
                                                   wfn_, options_);
      potential->compute();
      oepdev::ESPSolver esp(potential);
      esp.compute();
     
      for (int i=0; i<esp.charges()->nrow(); ++i) {
           oepTypes_["V"].matrix->set(i, 0, esp.charges()->get(i, 0));
      }
      psi::timer_off("OEPDEV: Electrostatic Energy OEP -> fitting ESP charges");

  } else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
  }

}
double ElectrostaticEnergyOEPotential::compute_3D_V(const double& x, const double& y, const double& z){
      double val = 0.0;
      // ===> Nuclear contribution <=== //
      for (int i=0; i<wfn_->molecule()->natom(); ++i) {
           val+= (double)wfn_->molecule()->Z(i) /
                     sqrt( pow(wfn_->molecule()->x(i)-x, 2.0) 
                         + pow(wfn_->molecule()->y(i)-y, 2.0) 
                         + pow(wfn_->molecule()->z(i)-z, 2.0) );
      }

      // ===> Electronic contribution <=== //
      double v;
      potInt_->set_charge_field(x, y, z);
      OEInt_ = potInt_;
      OEInt_->compute(potMat_);
      for (int i=0; i<primary_->nbf(); ++i) {
           for (int j=0; j<=i; ++j) {
                v= (wfn_->Da()->get(i,j) + wfn_->Db()->get(i,j)) * potMat_->get(i,j);
                val += 2.0*v;
                if (i==j) val -= v;
           }
      }
      potMat_->zero();
 return val;
}
void ElectrostaticEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) 
{
   double val;
   if (oepType == "V" || oepType == "TOTAL") val = compute_3D_V(x, y, z);
   else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
   }

   // Assign final value
   v->set(0, val);
}
void ElectrostaticEnergyOEPotential::print_header(void) const 
{
   psi::outfile->Printf("  ==> OEPotential: %s <==\n\n", name_.c_str());
   oepTypes_.at("V").matrix->print();
}


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

RepulsionEnergyOEPotential::~RepulsionEnergyOEPotential() {}
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
   //oepTypes_.push_back("Murrell-etal.S1");
   //oepTypes_.push_back("Otto-Ladik.S2");
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.S1") compute_murrell_etal_s1();
  else if (oepType ==   "Otto-Ladik.S2") compute_murrell_etal_s2();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void RepulsionEnergyOEPotential::compute_murrell_etal_s1() 
{
   // ---> Determine the target basis set for generalized density fitting <--- //
   std::shared_ptr<psi::BasisSet> target = intermediate_;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") target = auxiliary_;

   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> Ca_occ = wfn_->Ca_subset("AO","OCC");

   psi::IntegralFactory fact_1(primary_, target);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, target);
  
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());

   // ===> Compute One-Electron Integrals <=== //
   potInt->compute(Vao);

   // ===> Compute Nuclear Contribution to DF Vector <=== //
   std::shared_ptr<psi::Matrix> V = psi::Matrix::doublet(Vao, Ca_occ, true, false);

   // ===> Compute Electronic Contribution to DF Vector <=== //
   if (true) {
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
   }

   // ===> Perform The Generalized Density Fitting <=== // 
   std::shared_ptr<oepdev::GeneralizedDensityFit> gdf;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, V);}
   else                                                {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, intermediate_, V);}
   std::shared_ptr<psi::Matrix> G = gdf->compute();
   
   // ===> Save and Finish <=== //
   oepTypes_.at("Murrell-etal.S1").matrix->copy(G);
   //G->print();
}
void RepulsionEnergyOEPotential::compute_murrell_etal_s2() 
{
// TODO
}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void RepulsionEnergyOEPotential::print_header(void) const {}

// <============== CT Energy ===============> //
ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   throw psi::PSIEXCEPTION("OEPDEV: Construction of Charge Transfer Energy OEP requires auxiliary basis set!\n");
   common_init();
}

ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, 
                                                                 SharedBasisSet auxiliary, 
                                                                 SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

ChargeTransferEnergyOEPotential::~ChargeTransferEnergyOEPotential() {}
void ChargeTransferEnergyOEPotential::common_init() 
{
    int n1 = primary_->nbf();
    int n2 = auxiliary_->nbf();
    int n3 = wfn_->molecule()->natom();

    SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G ", n2, n1);
    SharedMatrix mat_2 = std::make_shared<psi::Matrix>("Q1", n3, n1);
    SharedMatrix mat_3 = std::make_shared<psi::Matrix>("Q2", n3, 1);

    OEPType type_1 = {"Otto-Ladik.V1", true , n1, mat_1};
    OEPType type_2 = {"Otto-Ladik.V2", false, n1, mat_2};
    OEPType type_3 = {"Otto-Ladik.V3", false,  1, mat_3};

    oepTypes_[type_1.name] = type_1;
    oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
}

void ChargeTransferEnergyOEPotential::compute(const std::string& oepType) {}
void ChargeTransferEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void ChargeTransferEnergyOEPotential::print_header(void) const {}

void ChargeTransferEnergyOEPotential::compute_otto_ladik_v1() {}
void ChargeTransferEnergyOEPotential::compute_otto_ladik_v2() {}
void ChargeTransferEnergyOEPotential::compute_otto_ladik_v3() {}

// <============== EET Coupling ==============> //

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

EETCouplingOEPotential::EETCouplingOEPotential(SharedWavefunction wfn, 
                                               SharedBasisSet auxiliary, 
                                               SharedBasisSet intermediate, Options& options) 
 : OEPotential(wfn, auxiliary, intermediate, options)
{ 
   common_init();
}

EETCouplingOEPotential::~EETCouplingOEPotential() {}
void EETCouplingOEPotential::common_init() 
{
    int n1 = primary_->nbf();
    int n2 = auxiliary_->nbf();
    int n3 = wfn_->molecule()->natom();

    // Implement these matrices (here provide correct dimensions)
    SharedMatrix mat_1 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_2 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  
    SharedMatrix mat_3 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_4 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  
    SharedMatrix mat_5 = std::make_shared<psi::Matrix>("ET1-NotImplemented", 1, 1);  
    SharedMatrix mat_6 = std::make_shared<psi::Matrix>("ET2-NotImplemented", 1, 1);  

    // Provide correct classes (DF- or ESP-based) and number of OEP's involved in each type
    OEPType type_1 = {"Fujimoto.ET1", false,  0, mat_1};
    OEPType type_2 = {"Fujimoto.ET2", false,  0, mat_2};
    OEPType type_3 = {"Fujimoto.HT1", false,  0, mat_3};
    OEPType type_4 = {"Fujimoto.HT2", false,  0, mat_4};
    OEPType type_5 = {"Fujimoto.CT1", false,  0, mat_5};
    OEPType type_6 = {"Fujimoto.CT2", false,  0, mat_6};

    oepTypes_[type_1.name] = type_1;
    oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
    oepTypes_[type_4.name] = type_4;
    oepTypes_[type_5.name] = type_5;
    oepTypes_[type_6.name] = type_6;
}

void EETCouplingOEPotential::compute(const std::string& oepType) {}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void EETCouplingOEPotential::print_header(void) const {}
