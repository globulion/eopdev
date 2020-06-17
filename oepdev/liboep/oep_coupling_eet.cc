#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;

// <============== EET Coupling ==============> //
SharedOEPotential EETCouplingOEPotential::clone(void) const {
    SharedOEPotential temp = std::make_shared<EETCouplingOEPotential>(this);
    return temp;
}
EETCouplingOEPotential::EETCouplingOEPotential(const EETCouplingOEPotential* f) 
 : OEPotential(f) {} 
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
    SharedMatrix mat_exch= std::make_shared<psi::Matrix>("G(EXCH)", n1, n1);
    SharedMatrix mat_ct_m= std::make_shared<psi::Matrix>("R(HL)", 1, 1);

    // Provide correct classes (DF- or ESP-based) and number of OEP's involved in each type
    OEPType type_1 = {"Fujimoto.GDF", true, 4, mat_gdf, nullptr, nullptr};
    OEPType type_2 = {"Fujimoto.CIS", false,1, std::make_shared<psi::Matrix>(), nullptr, nullptr};
    OEPType type_3 = {"Fujimoto.EXCH", false, n1, mat_exch, nullptr, nullptr};
    OEPType type_4 = {"Fujimoto.CT_M", false, 1, mat_ct_m, nullptr, nullptr};

    // Register OEPType's
    oepTypes_[type_1.name] = type_1;
    oepTypes_[type_2.name] = type_2;
    oepTypes_[type_3.name] = type_3;
    oepTypes_[type_4.name] = type_4;
}

void EETCouplingOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Fujimoto.GDF"   ) compute_fujimoto_gdf();
  else if (oepType == "Fujimoto.CIS"   ) compute_fujimoto_cis();
  else if (oepType == "Fujimoto.EXCH"  ) compute_fujimoto_exch();
  else if (oepType == "Fujimoto.CT_M"  ) compute_fujimoto_ct_m();
  else throw psi::PSIEXCEPTION("OEPDEV: Error. Incorrect OEP type specified!\n");
}
void EETCouplingOEPotential::compute_fujimoto_cis()
{
      const bool symm = options_.get_bool("TrCAMM_SYMMETRIZE");
      int I = options_.get_int("EXCITED_STATE");
      const int nH = options_.get_int("OEPDEV_SOLVER_EET_HOMO");
      const int nL = options_.get_int("OEPDEV_SOLVER_EET_LUMO");

      std::shared_ptr<CISComputer> cis = CISComputer::build("RESTRICTED", wfn_, options_, "RHF"); 
      cis->compute();
      cis->determine_electronic_state(I);
      oepTypes_.at("Fujimoto.CIS").cis_data = cis->data(I, nH, nL, symm);
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
   const int nH = options_.get_int("OEPDEV_SOLVER_EET_HOMO");
   const int nL = options_.get_int("OEPDEV_SOLVER_EET_LUMO");
   const int homo = wfn_->nalpha() - 1 - nH;
   const int lumo = nL;
   psi::SharedVector C_H = cOcc_->get_column(0, homo);
   psi::SharedVector C_L = cVir_->get_column(0, lumo);

   psi::IntegralFactory fact_1(primary_, target, primary_, target);
   psi::IntegralFactory fact_2(primary_, primary_, primary_, target);

   std::shared_ptr<psi::OneBodyAOInt> kinInt(fact_1.ao_kinetic());
   std::shared_ptr<psi::OneBodyAOInt> potInt(fact_1.ao_potential());

   // ===> Compute One-Electron Integrals <=== //
   kinInt->compute(Tao); Tao->scale(0.5);
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
             int i = intsIter->i();  // \gamma     : n
             int j = intsIter->j();  // \delta     : n
             int k = intsIter->k();  // \beta      : n
             int l = intsIter->l();  // \alpha     : Q

             double eri = buffer[intsIter->index()];

             double dij = d[i][j]; double dik = d[i][k];
             double cli = cl[i]; double clk = cl[k]; double clj = cl[j];
             double chi = ch[i]; double chk = ch[k]; double chj = ch[j];

             double vl_l = eri * (2.0*clk*dij - clj*dik);
             double vl_h =-eri * (2.0*chk*dij - chj*dik);      // OK
             double vl_et= eri * (2.0*chk*cli - clk*chi)* chj;
             double vl_ht= eri * (2.0*clk*chi - chk*cli)* clj; // OK

             // ET (L)
             V[l][0] += vl_l;
             // ET (HL)
             V[l][1] += vl_l + vl_et;
             // HT (H)
             V[l][2] += vl_h;
             // HT (HL)
             V[l][3] += vl_h + vl_ht;
        }
   }

   // ===> Perform The Generalized Density Fitting <=== // 
   std::shared_ptr<oepdev::GeneralizedDensityFit> gdf;
   if (options_.get_str("OEPDEV_DF_TYPE") == "SINGLE") {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, vao);}
   else                                                {gdf = oepdev::GeneralizedDensityFit::build(auxiliary_, intermediate_,vao);}
   std::shared_ptr<psi::Matrix> G = gdf->compute();
   
   // ===> Save and Finish <=== //
   oepTypes_.at("Fujimoto.GDF").matrix = G;
   if (options_.get_int("PRINT") > 1) G->print();

}
void EETCouplingOEPotential::compute_fujimoto_exch()
{
   // ===> Allocate <=== //
   std::shared_ptr<psi::Matrix> G = std::make_shared<psi::Matrix>("G" , primary_->nbf(), primary_->nbf());
   psi::IntegralFactory fact(primary_, primary_, primary_, primary_);

   // Compute //
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact.eri());
   const double * buffer = tei->buffer();

   std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact, "ALL");

   double** g = G->pointer(); 
   for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
   {
        shellIter->compute_shell(tei);
        std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
        for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
        {
             int i = intsIter->i();
             int j = intsIter->j();
             int k = intsIter->k();
             int l = intsIter->l();

             double eri = buffer[intsIter->index()];
             if ((i==j)&&(k==l)) g[i][j] = eri;
        }
   }

   // ===> Save and Finish <=== //
   oepTypes_.at("Fujimoto.EXCH").matrix = G;
   if (options_.get_int("PRINT") > 1) G->print();
}
void EETCouplingOEPotential::compute_fujimoto_ct_m()
{
   // ===> Allocate <=== //
   psi::IntegralFactory fact(primary_, primary_, primary_, primary_);

   // Compute //
   std::shared_ptr<psi::TwoBodyAOInt> tei(fact.eri());
   const double * buffer = tei->buffer();

   std::shared_ptr<oepdev::ShellCombinationsIterator> shellIter = oepdev::ShellCombinationsIterator::build(fact, "ALL");

   double v = 0.0; // (HH|LL) integral
   const int nH = options_.get_int("OEPDEV_SOLVER_EET_HOMO");
   const int nL = options_.get_int("OEPDEV_SOLVER_EET_LUMO");
   const int homo = wfn_->nalpha() - 1 - nH;
   const int lumo = nL;
   psi::SharedVector Vh = cOcc_->get_column(0, homo);
   psi::SharedVector Vl = cVir_->get_column(0, lumo);
   double* ch = Vh->pointer();
   double* cl = Vl->pointer();
 //double* ch = cOcc_->get_column(0, homo)->pointer(); --> memory leak!!!
 //double* cl = cVir_->get_column(0, lumo)->pointer(); --> memory leak!!!
   for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
   {
        shellIter->compute_shell(tei);
        std::shared_ptr<oepdev::AOIntegralsIterator> intsIter = shellIter->ao_iterator("ALL");
        for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
        {
             int i = intsIter->i();
             int j = intsIter->j();
             int k = intsIter->k();
             int l = intsIter->l();

             v += ch[i] * ch[j] * cl[k] * cl[l] * buffer[intsIter->index()];
        }
   }

   // ===> Save and Finish <=== //
   oepTypes_.at("Fujimoto.CT_M").matrix->set(0, 0, v);
   if (options_.get_int("PRINT") > 1) psi::outfile->Printf("\n (HH|LL) = %14.6f\n\n", v);
}

void EETCouplingOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) {}
void EETCouplingOEPotential::print_header(void) const {}

void EETCouplingOEPotential::rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) {

  // Potential "Fujimoto.GDF"
  psi::SharedMatrix Rj = R_aux->clone(); Rj->invert(); Rj->transpose_this();
  psi::SharedMatrix new_matrix = psi::Matrix::doublet(Rj, oepTypes_.at("Fujimoto.GDF").matrix, true, false);
  oepTypes_.at("Fujimoto.GDF").matrix = new_matrix;

  // Potential "Fujimoto.EXCH" -> TODO
  // add this to rotate oepTypes_.at("Fujimoto.CIS").matrix !

  // Potential "Fujimoto.CIS"
  psi::SharedMatrix Ri = R_prim->clone(); Ri->invert(); Ri->transpose_this();
  psi::SharedMatrix new_Pe = psi::Matrix::triplet(Ri, oepTypes_.at("Fujimoto.CIS").cis_data->Pe , Ri, true, false, false);
  psi::SharedMatrix new_Peg= psi::Matrix::triplet(Ri, oepTypes_.at("Fujimoto.CIS").cis_data->Peg, Ri, true, false, false);
  oepTypes_.at("Fujimoto.CIS").cis_data->Pe = new_Pe;
  oepTypes_.at("Fujimoto.CIS").cis_data->Peg = new_Peg;
  oepTypes_.at("Fujimoto.CIS").cis_data->trcamm->rotate(r);
  oepTypes_.at("Fujimoto.CIS").cis_data->camm_homo->rotate(r);
  oepTypes_.at("Fujimoto.CIS").cis_data->camm_lumo->rotate(r);

  this->rotate_basic(r, R_prim, R_aux);
}
void EETCouplingOEPotential::translate_oep(psi::SharedVector t) {
  // Potential "Fujimoto.CIS"
  oepTypes_.at("Fujimoto.CIS").cis_data->trcamm->translate(t);
  oepTypes_.at("Fujimoto.CIS").cis_data->camm_homo->translate(t);
  oepTypes_.at("Fujimoto.CIS").cis_data->camm_lumo->translate(t);

  this->translate_basic(t);
}
