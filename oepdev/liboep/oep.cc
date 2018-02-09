#include "oep.h"
#include "../libutil/space3d.h"
#include "../libutil/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"


using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::ScalarField3D>;

OEPotential::OEPotential(SharedWavefunction wfn, Options& options) 
   : wfn_(wfn),
     options_(options),
     is_density_fitted(false),
     is_esp_based(true)
{
    common_init();
}

OEPotential::OEPotential(SharedWavefunction wfn, SharedBasisSet auxiliary, Options& options) 
   : wfn_(wfn),
     auxiliary_(auxiliary),
     options_(options),
     is_density_fitted(true),
     is_esp_based(false)
{
    common_init();
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
                                                                           SharedBasisSet auxiliary, Options& options)
{
   std::shared_ptr<OEPotential> oep;

   if      (category == "ELECTROSTATIC ENERGY"  )  throw PSIEXCEPTION("OEPDEV: OEPotential build. DF-based OEP not available for Electrostatic Energy!");
   else if (category == "REPULSION ENERGY"      )  oep = std::make_shared< RepulsionEnergyOEPotential>(wfn, auxiliary, options);
   else if (category == "CHARGE TRANSFER ENERGY")  oep = std::make_shared< ChargeTransferEnergyOEPotential>(wfn, auxiliary, options);
   else if (category == "EET COUPLING CONSTANT" )  oep = std::make_shared<    EETCouplingOEPotential>(wfn, auxiliary, options);
   else  
            throw PSIEXCEPTION("OEPDEV: OEPotential build. Unrecognized OEP category.");

   return oep;
}



OEPotential::~OEPotential() {}

void OEPotential::common_init(void) 
{
   name_        = "default";
   primary_     = wfn_->basisset();
   intsFactory_ = std::make_shared<psi::IntegralFactory>(primary_, primary_);
   potMat_      = std::make_shared<psi::Matrix>("Potential Integrals", primary_->nbf(), primary_->nbf());
   potInt_      = std::make_shared<oepdev::PotentialInt>(intsFactory_->spherical_transform(), primary_, primary_, 0);
}

void OEPotential::compute(const std::string& oepType) {}
void OEPotential::compute(void) { for ( const std::string& type : oepTypes_ ) this->compute(type); }
void OEPotential::write_cube(const std::string& oepType, const std::string& fileName) 
{
   OEPotential3D<OEPotential> cube(60, 60, 60, 10.0, 10.0, 10.0, shared_from_this(), oepType, options_);
   cube.write_cube_file(fileName);
}
void OEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
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
   oepTypes_.push_back("V");
}

void ElectrostaticEnergyOEPotential::compute(const std::string& oepType) 
{
  if (oepType == "V" || oepType == "TOTAL") {

      psi::timer_on("OEPDEV: Electrostatic Energy OEP -> fitting ESP charges");
      SharedField3D potential = oepdev::ScalarField3D::build("ELECTROSTATIC POTENTIAL", 
                                                   options_.get_int   ("ESP_NPOINTS_PER_ATOM") * wfn_->molecule()->natom(), 
                                                   options_.get_double("ESP_PAD_SPHERE"      ), 
                                                   wfn_, options_);
      potential->compute();
      oepdev::ESPSolver esp(potential);
      esp.compute();
     
      oepMatrices_["V"] = std::make_shared<Matrix>("V", esp.charges()->dim(), 1);
      for (int i=0; i<esp.charges()->dim(); ++i) {
           oepMatrices_["V"]->set(i, 0, esp.charges()->get(i));
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
void ElectrostaticEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) 
{
   double val;
   if (oepType == "V" || oepType == "TOTAL") val = compute_3D_V(x, y, z);
   else {
      throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
   }

   // Assign final value
   v = val;
}
void ElectrostaticEnergyOEPotential::print_header(void) const 
{
   psi::outfile->Printf("  ==> OEPotential: %s <==\n\n", name_.c_str());
   oepMatrices_.at("V")->print();
}


// <============== Repulsion Energy ==============> //

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(SharedWavefunction wfn, 
                                                       SharedBasisSet auxiliary, Options& options) 
 : OEPotential(wfn, auxiliary, options)
{ 
   common_init();
}

RepulsionEnergyOEPotential::~RepulsionEnergyOEPotential() {}
void RepulsionEnergyOEPotential::common_init() 
{
   oepTypes_.push_back("Murrell-etal.S1");
   oepTypes_.push_back("Otto-Ladik.S2");
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.S1") compute_murrell_etal_s1();
  else if (oepType ==   "Otto-Ladik.S2") compute_murrell_etal_s2();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void RepulsionEnergyOEPotential::compute_murrell_etal_s1() 
{
   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao" , primary_->nbf(), auxiliary_->nbf());
   std::shared_ptr<psi::Matrix> Sao = std::make_shared<psi::Matrix>("Sao" , auxiliary_->nbf(), auxiliary_->nbf());
   std::shared_ptr<psi::Matrix> Ca_occ = wfn_->Ca_subset("AO","OCC");

   psi::IntegralFactory fact_a(auxiliary_, auxiliary_);
   psi::IntegralFactory fact_1(primary_, auxiliary_);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, auxiliary_);
  
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());
   std::shared_ptr<psi::OneBodyAOInt> ovlInt(fact_a.ao_overlap());

   // ===> Compute One-Electron Integrals <=== //
   potInt->compute(Vao);
   ovlInt->compute(Sao);

   //psi::shared_ptr<psi::Matrix> I = psi::Matrix::doublet(Sao, Sao, false, false);
   //I->print();

   // ===> Compute Nuclear Contribution to DF Vector <=== //
   std::shared_ptr<psi::Matrix> V = psi::Matrix::doublet(Vao, Ca_occ, true, false);

   // ===> Compute Electronic Contribution to DF Vector <=== //
   if (true) {
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact_2.eri());                                           
   const double * buffer = tei->buffer();
                                                                                                  
   oepdev::AllAOShellCombinationsIterator shellIter(fact_2);
   int i, j, k, l;
   double integral, dij, dik;
   
   double** v = V->pointer();  
   double** c = Ca_occ->pointer();
   double** d = wfn_->Da()->pointer();                                                                                 
   for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
   {
        shellIter.compute_shell(tei);
        oepdev::AllAOIntegralsIterator intsIter(shellIter);
        for (intsIter.first(); intsIter.is_done() == false; intsIter.next())
        {
             i = intsIter.i();  // \mu      : n
             j = intsIter.j();  // \nu      : n
             k = intsIter.k();  // \alpha   : n
             l = intsIter.l();  // \xi      : Q

             integral = buffer[intsIter.index()];

             dij = d[i][j]; dik = d[i][k];
             for (int a = 0; a < wfn_->doccpi()[0]; ++a) 
             {
                v[l][a] += (2.0 * c[k][a]*dij - c[j][a]*dik) * integral;
             }
        }
   }
   }

   //// ===> Perform Generalized Density Fitting <=== // 
   Sao->invert();
   std::shared_ptr<psi::Matrix> G = psi::Matrix::doublet(Sao, V, false, false);
   G->set_name("G(S^{-1})");
   //double** pS = Sao->pointer();
   //double** v = V->pointer();  
   //int Q = auxiliary_->nbf();
   //int N = wfn_->doccpi()[0];
   //int* ipiv = init_int_array(Q);
   ////C_DGESV(Q, N, &(pS[0][0]), Q, &(ipiv[0]), &(v[0][0]), Q);
   //int info = C_DGESV(Q, N, Sao->pointer()[0], Q, ipiv, V->pointer()[0], Q);
   //cout << info << endl;

   //free(ipiv);
   //V->print();
   
   // ===> Save and Finish <=== //
   oepMatrices_["Murrell-etal.S1"] = G;
   //G->print();
}
void RepulsionEnergyOEPotential::compute_murrell_etal_s2() 
{
}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
void RepulsionEnergyOEPotential::print_header(void) const {}

// <============== CT Energy ===============> //
ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, Options& options) 
 : OEPotential(wfn, options)
{ 
   throw psi::PSIEXCEPTION("OEPDEV: Construction of Charge Transfer Energy OEP requires auxiliary basis set!\n");
   common_init();
}

ChargeTransferEnergyOEPotential::ChargeTransferEnergyOEPotential(SharedWavefunction wfn, 
                                                                 SharedBasisSet auxiliary, Options& options) 
 : OEPotential(wfn, auxiliary, options)
{ 
   common_init();
}

ChargeTransferEnergyOEPotential::~ChargeTransferEnergyOEPotential() {}
void ChargeTransferEnergyOEPotential::common_init() 
{
    oepTypes_.push_back("Otto-Ladik.V1");
    oepTypes_.push_back("Otto-Ladik.V2");
    oepTypes_.push_back("Otto-Ladik.V3");
}

void ChargeTransferEnergyOEPotential::compute(const std::string& oepType) {}
void ChargeTransferEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
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
                                               SharedBasisSet auxiliary, Options& options) 
 : OEPotential(wfn, auxiliary, options)
{ 
   common_init();
}

EETCouplingOEPotential::~EETCouplingOEPotential() {}
void EETCouplingOEPotential::common_init() 
{
    oepTypes_.push_back("Fujimoto.ET1");
    oepTypes_.push_back("Fujimoto.ET2");
    oepTypes_.push_back("Fujimoto.HT1");
    oepTypes_.push_back("Fujimoto.HT2");
    oepTypes_.push_back("Fujimoto.CT1");
    oepTypes_.push_back("Fujimoto.CT2");
}

void EETCouplingOEPotential::compute(const std::string& oepType) {}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, double& v) {}
void EETCouplingOEPotential::print_header(void) const {}
