#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../lib3d/dmtp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;
using SharedDMTPole = std::shared_ptr<oepdev::DMTPole>;

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
    name_ = "HF Charge-Transfer Energy";

    int n1 = primary_->nbf();
    int n2 = auxiliary_->nbf();
    int n3 = wfn_->molecule()->natom();

    SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G ", n2, n1);
  //SharedMatrix mat_2 = std::make_shared<psi::Matrix>("Q1", n3, n1);
    SharedMatrix mat_3 = std::make_shared<psi::Matrix>("Q2", n3, 1);

    OEPType type_1 = {"Murrell-etal.V1.GDF", true , n1, mat_1};
  //OEPType type_2 = {"Murrell-etal.V2", false, n1, mat_2};
    OEPType type_3 = {"Murrell-etal.V3.CAMM-nj", false, n1, std::make_shared<psi::Matrix>()};

    oepTypes_[type_1.name] = type_1;
  //oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
}

void ChargeTransferEnergyOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.V1.GDF"    ) compute_murrell_etal_v1_gdf    ();
  else if (oepType == "Murrell-etal.V3.CAMM-nj") compute_murrell_etal_v3_camm_nj();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void ChargeTransferEnergyOEPotential::compute_murrell_etal_v1_gdf()
{
   // ---> Determine the target basis set for generalized density fitting <--- //
   std::shared_ptr<psi::BasisSet> target = intermediate_;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") target = auxiliary_;

   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> Ca_vir = wfn_->Ca_subset("AO","VIR");
   const int nvir = Ca_vir->ncol();

   psi::IntegralFactory fact_1(primary_, target, primary_, target);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, target);
  
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());

   // ===> Compute One-Electron Integrals <=== //
   potInt->compute(Vao);

   // ===> Compute Nuclear Contribution to DF Vector <=== //
   std::shared_ptr<psi::Matrix> V = psi::Matrix::doublet(Vao, Ca_vir, true, false);

   // ===> Compute Electronic Contribution to DF Vector <=== //
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact_2.eri());                                           
   const double * buffer = tei->buffer();

   std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact_2, "ALL");
   int i, j, k, l;
   double integral, dij, dik;
 
   double** v = V->pointer();  
   double** c = Ca_vir->pointer();
   double** d = wfn_->Da()->pointer();
   for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
   {
        shellIter->compute_shell(tei);
        std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
        for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
        {
             i = intsIter->i();  // \beta      : n
             j = intsIter->j();  // \gamma     : n
             k = intsIter->k();  // \alpha     : n
             l = intsIter->l();  // \xi        : Q

             integral = buffer[intsIter->index()];

             dij = d[i][j]; dik = d[i][k];
             for (int a = 0; a < nvir; ++a) 
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
   oepTypes_.at("Murrell-etal.V1.GDF").matrix->copy(G);
   if (options_.get_int("PRINT") > 1) G->print();
}
void ChargeTransferEnergyOEPotential::compute_murrell_etal_v3_camm_nj()
{ 
  // ==> Sizing <== //
  int nocc = cOcc_->ncol();
  int nvir = cVir_->ncol();
  int nbf = primary_->nbf();

  // ==> Initialize CAMM object <== //
  SharedDMTPole camm = oepdev::DMTPole::build("CAMM", wfn_, nocc*nvir);

  // TODO
}


//void ChargeTransferEnergyOEPotential::compute_murrell_etal_v2() {}

void ChargeTransferEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void ChargeTransferEnergyOEPotential::print_header(void) const {}


