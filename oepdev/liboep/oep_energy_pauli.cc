#include "oep.h"
#include "oep_gdf.h"
#include "../lib3d/esp.h"
#include "../lib3d/dmtp.h"
#include "../libutil/integrals_iter.h"
#include "psi4/libqt/qt.h"

using namespace oepdev;

using SharedField3D = std::shared_ptr<oepdev::Field3D>;
using SharedDMTPole = std::shared_ptr<oepdev::DMTPole>;

// <============== Repulsion Energy ==============> //
SharedOEPotential RepulsionEnergyOEPotential::clone(void) const {
    SharedOEPotential temp = std::make_shared<RepulsionEnergyOEPotential>(this);
    return temp;
}
RepulsionEnergyOEPotential::RepulsionEnergyOEPotential(const RepulsionEnergyOEPotential* f) 
 : OEPotential(f)  {vec_otto_ladik_s2_ = new double[1];}
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
  if (vec_otto_ladik_s2_) {
     delete[] vec_otto_ladik_s2_;
   }
}
void RepulsionEnergyOEPotential::common_init() 
{
   name_ = "HF Repulsion Energy";

   int n1 = wfn_->Ca_subset("AO","OCC")->ncol();
   int n2 = auxiliary_->nbf();
   int n3 = wfn_->molecule()->natom();

   SharedMatrix mat_1 = std::make_shared<psi::Matrix>("G(S^{-1})", n2, n1);
   SharedMatrix mat_2 = std::make_shared<psi::Matrix>("G(S^{-2})", n3, n1);

   OEPType type_1 = {"Murrell-etal.S1"       , true , n1, mat_1, nullptr, nullptr};
   OEPType type_2 = {"Otto-Ladik.S2.ESP"     , false, n1, mat_2, nullptr, nullptr};
   OEPType type_3 = {"Otto-Ladik.S2.CAMM.a"  , false, n1, std::make_shared<psi::Matrix>(), nullptr, nullptr};
   OEPType type_4 = {"Otto-Ladik.S2.CAMM.A"  , false, n1, std::make_shared<psi::Matrix>(), nullptr, nullptr};

   oepTypes_[type_1.name] = type_1; 
   oepTypes_[type_2.name] = type_2;
   oepTypes_[type_3.name] = type_3;
   oepTypes_[type_4.name] = type_4;

   //
   vec_otto_ladik_s2_ = new double[n1];
}

void RepulsionEnergyOEPotential::compute(const std::string& oepType) 
{
  if      (oepType == "Murrell-etal.S1"     ) compute_murrell_etal_s1()      ;
  else if (oepType == "Otto-Ladik.S2.ESP"   ) compute_otto_ladik_s2_esp()    ;
  else if (oepType == "Otto-Ladik.S2.CAMM.a") compute_otto_ladik_s2_camm_a() ;
  else if (oepType == "Otto-Ladik.S2.CAMM.A") compute_otto_ladik_s2_camm_A() ;
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
   std::shared_ptr<psi::Matrix> Ca_occ;
   //if (options_.get_bool("OEPDEV_OEP_REP_LOCALIZE")) {
   if (this->use_localized_orbitals) {
       if (!localizer_) this->localize();
       Ca_occ = this->lOcc_;
   }
   else {Ca_occ = this->cOcc_;} //wfn_->Ca_subset("AO","OCC");

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
   oepTypes_.at("Murrell-etal.S1").matrix = G;
   if (options_.get_int("PRINT") > 1) G->print();
   //psi::timer_off("OEP    E(Paul) Murrell-etal S1  ");
}
void RepulsionEnergyOEPotential::compute_otto_ladik_s2_esp() 
{
      //psi::timer_on("OEP    E(Paul) Otto-Ladik.S2.ESP    ");
      std::shared_ptr<OEPotential3D<OEPotential>> oeps3d = this->make_oeps3d("Otto-Ladik.S2.ESP");
      oeps3d->compute();
      ESPSolver esp(oeps3d);
      esp.set_charge_sums(0.5);
      esp.compute();
     
      for (int i=0; i<esp.charges()->nrow(); ++i) {
      for (int o=0; o<oepTypes_["Otto-Ladik.S2.ESP"].n; ++o) {
           oepTypes_["Otto-Ladik.S2.ESP"].matrix->set(i, o, esp.charges()->get(i, o));
      }}
      if (options_.get_int("PRINT") > 1) esp.charges()->print();
      //psi::timer_off("OEP    E(Paul) Otto-Ladik.S2.ESP    ");
}
void RepulsionEnergyOEPotential::compute_otto_ladik_s2_camm_a() 
{
   // ==> Sizing <== //
   int nocc = oepTypes_["Otto-Ladik.S2.CAMM.a"].n;
   int nbf = primary_->nbf();

   // ==> Initialize CAMM object <== //
   SharedDMTPole camm = oepdev::DMTPole::build("CAMM", wfn_, nocc);

   // ==> Choose which orbitals to use <== //
   psi::SharedMatrix Ca_occ;
   //if (options_.get_bool("OEPDEV_OEP_REP_LOCALIZE")) {
   if (this->use_localized_orbitals) {
       if (!localizer_) this->localize();
       Ca_occ = this->lOcc_;
   }
   else {Ca_occ = this->cOcc_;} 

   // ==> Compute the vector of OED's <== //
   std::vector<psi::SharedMatrix> oeds;
   std::vector<bool> trans;
   for (int i=0; i<nocc; ++i) {

        // ==> Exclude nuclear part <== //
	trans.push_back(true);

	// ==> Set up OED <== //
	psi::SharedMatrix oed = std::make_shared<psi::Matrix>("", nbf, nbf);
	double** oed_p = oed->pointer();

	for (int a=0; a<nbf; ++a) {
        for (int b=0; b<nbf; ++b) {
	     oed_p[a][b] = Ca_occ->get(a,i) * Ca_occ->get(b,i);
        }
	}
	oeds.push_back(oed);
   }

   // ==> Compute CAMM sets <== //
   camm->compute(oeds, trans);

   // ==> Save <== //
   oepTypes_["Otto-Ladik.S2.CAMM.a"].dmtp = camm;

}
void RepulsionEnergyOEPotential::compute_otto_ladik_s2_camm_A() 
{
   // ==> Sizing <== //
   int nocc = oepTypes_["Otto-Ladik.S2.CAMM.A"].n;
   int nbf = primary_->nbf();

   // ==> Initialize CAMM object <== //
   SharedDMTPole camm = oepdev::DMTPole::build("CAMM", wfn_, nocc);
   SharedMatrix Da = wfn_->Da()->clone();

   // ==> Choose which orbitals to use <== //
   psi::SharedMatrix Ca_occ;
   //if (options_.get_bool("OEPDEV_OEP_REP_LOCALIZE")) {
   if (this->use_localized_orbitals) {
       if (!localizer_) this->localize();
       Ca_occ = this->lOcc_;
   }
   else {Ca_occ = this->cOcc_;} 

   // ==> Compute the vector of OED's <== //
   std::vector<psi::SharedMatrix> oeds;
   std::vector<bool> trans;
   for (int i=0; i<nocc; ++i) {

        // ==> Include nuclear part <== //
	trans.push_back(false);

	// ==> Set up OED <== //
	psi::SharedMatrix oed = std::make_shared<psi::Matrix>("", nbf, nbf);
	double** oed_p = oed->pointer();

	for (int a=0; a<nbf; ++a) {
        for (int b=0; b<nbf; ++b) {
	     oed_p[a][b] = Ca_occ->get(a,i) * Ca_occ->get(b,i);
        }
	}
	oed->scale(-0.5);
	oed->axpy(2.0, Da);
	oeds.push_back(oed);
   }

   // ==> Compute CAMM sets <== //
   camm->compute(oeds, trans);

   // ==> Save <== //
   oepTypes_["Otto-Ladik.S2.CAMM.A"].dmtp = camm;
}
void RepulsionEnergyOEPotential::compute_3D(const std::string& oepType, const double& x, const double& y, const double& z, std::shared_ptr<psi::Vector>& v) 
{
   if (oepType == "Otto-Ladik.S2.ESP") {
       this->compute_3D_otto_ladik_s2(x, y, z);
       // Assign final value
       for (int o = 0; o < oepTypes_["Otto-Ladik.S2.ESP"].n; ++o) v->set(o, vec_otto_ladik_s2_[o]);
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

  // ==> Choose which orbitals to use <== //
  psi::SharedMatrix Ca_occ;
  //if (options_.get_bool("OEPDEV_OEP_REP_LOCALIZE")) {
  if (this->use_localized_orbitals) {
      if (!localizer_) this->localize();
      Ca_occ = this->lOcc_;
  }
  else {Ca_occ = this->cOcc_;} 

  psi::SharedMatrix potMO = psi::Matrix::triplet(Ca_occ, potMat_, Ca_occ, true, false, false);
  potMat_->zero();
  //for (int o=0; o<oepTypes_["Otto-Ladik.S2.ESP"].n; ++o) val += 2.0 * potMO->get(o, o);
  val += 2.0 * potMO->trace();

  for (int o=0; o<oepTypes_["Otto-Ladik.S2.ESP"].n; ++o) {
     vec_otto_ladik_s2_[o] = val - 0.5 * potMO->get(o, o);
  }
}
void RepulsionEnergyOEPotential::print_header(void) const 
{
	psi::outfile->Printf("\n ===> Repulsion Energy OEPotential <=== \n");
	psi::outfile->Printf(  " ===>           HF level           <=== \n\n");
        psi::outfile->Printf(  "      S-1 term: Murrell et.al\n");
        psi::outfile->Printf(  "      S-2 term: Otto and Ladik\n");
}
void RepulsionEnergyOEPotential::rotate_oep(psi::SharedMatrix r, psi::SharedMatrix R_prim, psi::SharedMatrix R_aux) {

  // Potential "Murrell-etal.S1"
  psi::SharedMatrix Ri = R_aux->clone(); Ri->invert(); Ri->transpose_this();
  psi::SharedMatrix new_matrix = psi::Matrix::doublet(Ri, oepTypes_.at("Murrell-etal.S1").matrix, true, false);
  oepTypes_.at("Murrell-etal.S1").matrix = new_matrix;
//oepTypes_.at("Murrell-etal.S1").matrix->copy(new_matrix); --> wrong!!! why???
//oepTypes_["Murrell-etal.S1"].matrix->copy(new_matrix.get()); ---> wrong!!! why???

  // Potential "Otto-Ladik.S2.CAMM.a"
  oepTypes_.at("Otto-Ladik.S2.CAMM.a").dmtp->rotate(r);

  // Potential "Otto-Ladik.S2.CAMM.A"
  oepTypes_.at("Otto-Ladik.S2.CAMM.A").dmtp->rotate(r);

  this->rotate_basic(r, R_prim, R_aux);
}
void RepulsionEnergyOEPotential::translate_oep(psi::SharedVector t) {

  // Potential "Otto-Ladik.S2.CAMM.a"
  oepTypes_.at("Otto-Ladik.S2.CAMM.a").dmtp->translate(t);
  // Potential "Otto-Ladik.S2.CAMM.A"
  oepTypes_.at("Otto-Ladik.S2.CAMM.A").dmtp->translate(t);

  this->translate_basic(t);
}
