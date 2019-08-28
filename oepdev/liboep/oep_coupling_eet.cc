#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

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

    // OEP Matrices
    SharedMatrix mat_gdf = std::make_shared<psi::Matrix>("ET and HT", n2, 4);  

    // Provide correct classes (DF- or ESP-based) and number of OEP's involved in each type
    OEPType type_1 = {"Fujimoto.GDF", true, 4, mat_gdf};

    // Register OEPType's
    oepTypes_[type_1.name] = type_1;
}

void EETCouplingOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Fujimoto.GDF"   ) compute_fujimoto_gdf();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void EETCouplingOEPotential::compute_fujimoto_gdf()
{
   // ---> Determine the target basis set for generalized density fitting <--- //
   std::shared_ptr<psi::BasisSet> target = intermediate_;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") target = auxiliary_;

   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> Tao = std::make_shared<psi::Matrix>("Tao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> Vao = std::make_shared<psi::Matrix>("Vao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> Hao = std::make_shared<psi::Matrix>("Hao" , primary_->nbf(), target->nbf());
   std::shared_ptr<psi::Matrix> vao = std::make_shared<psi::Matrix>("vao" , target->nbf(), 4);
   const int homo = wfn_->nalpha() - 1;
   const int lumo = 0;
   psi::SharedVector C_H = cOcc_->get_column(0, homo);
   psi::SharedVector C_L = cVir_->get_column(0, lumo);

   psi::IntegralFactory fact_1(primary_, target, primary_, target);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, target);

   std::shared_ptr<psi::OneBodyAOInt> kinInt(fact_1.ao_kinetic());
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());

   // ===> Compute One-Electron Integrals <=== //
   kinInt->compute(Tao);
   potInt->compute(Vao);

   // ===> Compute Hcore Matrix in Target/Primary Space <=== //
   Hao->add(Tao); Tao.reset();
   Hao->add(Vao); Vao.reset();
  
   // ===> Compute Hcore Matrix Contribution to DF Vectors <=== //
   double** V = vao->pointer();
   double** H = Hao->pointer();
   double* ch = C_H->pointer();
   double* cl = C_L->pointer();
   for (int i=0; i<target->nbf(); ++i) {
        double vl = 0.0;
        double vh = 0.0;
        for (int j=0; j<primary_->nbf(); ++j) {
             vl += cl[j] * H[j][i];
             vh += ch[j] * H[j][i];
        }
        for (int z=0; z<2; ++z) V[i][z] = vl;
        for (int z=2; z<4; ++z) V[i][z] =-vh;
   }
   Hao.reset();

   // ===> Compute Electronic Contribution to DF Vectors <=== //
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact_2.eri());                                           
   const double * buffer = tei->buffer();

   std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact_2, "ALL");
 
   double** d = wfn_->Da()->pointer();
   for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
   {
        shellIter->compute_shell(tei);
        std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
        for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
        {
             int i = intsIter->i();  // \beta      : n
             int j = intsIter->j();  // \gamma     : n
             int k = intsIter->k();  // \delta     : n
             int l = intsIter->l();  // \xi        : Q

             double eri = buffer[intsIter->index()];

             double dkj = d[k][j]; double dij = d[i][j];
             double cli = cl[i]; double clk = cl[k]; double clj = cl[j];
             double chi = ch[i]; double chk = ch[k]; double chj = ch[j];

             double vl_l = eri * (2.0*cli*dkj - clk*dij);
             double vl_h =-eri * (2.0*chi*dkj - chk*dij);
             double vl_et= eri * (2.0*chi * clj - cli * chj) * chk;
             double vl_ht= eri * (2.0*cli * chj - chi * clj) * clk;

             // ET A
             V[l][0] += vl_l + vl_et;
             // ET B
             V[l][1] += vl_l;
             // HT A
             V[l][2] += vl_h + vl_ht;
             // HT B
             V[l][3] += vl_h;
        }
   }

   // ===> Perform The Generalized Density Fitting <=== // 
   std::shared_ptr<oepdev::GeneralizedDensityFit> gdf;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, vao);}
   else                                                {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, intermediate_,vao);}
   std::shared_ptr<psi::Matrix> G = gdf->compute();
   
   // ===> Save and Finish <=== //
   oepTypes_.at("Fujimoto.GDF").matrix->copy(G);
   if (options_.get_int("PRINT") > 1) G->print();

}
void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void EETCouplingOEPotential::print_header(void) const {}
