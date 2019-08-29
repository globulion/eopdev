#include "ti_data.h"

using namespace std;
using namespace psi;
using namespace oepdev;

TIData::TIData()
 : s12(0.0), s13(0.0), s14(0.0), s23(0.0), s34(0.0), s24(0.0), de1(0.0), de2(0.0), v0({}),
   diagonal_correction(true), 
   mulliken_approximation(false),
   overlap_correction(true),
   trcamm_approximation(false),
   trcamm_convergence(oepdev::MultipoleConvergence::ConvergenceLevel::R5),
   c_(1.0)
{
}
TIData::~TIData() {}

void TIData::set_output_coupling_units_converter(double c) {this->c_ = c;}
void TIData::set_s(double s12, double s13, double s14, double s23, double s24, double s34)
{
 this->s12 = s12;
 this->s13 = s13;
 this->s14 = s14;
 this->s23 = s23;
 this->s24 = s24;
 this->s34 = s34;
}
void TIData::set_e(double e1, double e2, double e3, double e4)
{
 this->e1 = e1;
 this->e2 = e2;
 this->e3 = e3;
 this->e4 = e4;
}
void TIData::set_de(double de1, double de2)
{
 this->de1=de1;
 this->de2=de2;
}
double TIData::overlap_corrected(const std::string& type) {
  double v; bool direct = false; double s;
  if      (type == "COUL") {v = this->v0.at("COUL"); direct=true; s=this->s12;}
  else if (type == "EXCH") {v = this->v0.at("EXCH"); direct=true; s=this->s12;}
  else if (type == "OVRL") {v = 0.0; s= this->s12;} 
  else if (type == "ET1" ) {v = this->v0.at("ET1"); s=this->s13;}
  else if (type == "ET2" ) {v = this->v0.at("ET2"); s=this->s24;}
  else if (type == "HT1" ) {v = this->v0.at("HT1"); s=this->s14;}
  else if (type == "HT2" ) {v = this->v0.at("HT2"); s=this->s23;}
  else if (type == "CT"  ) {v = this->v0.at("CT" ); s=this->s34;}
  else if (type == "CT_M") {v = this->v0.at("CT_M");s=this->s34;}
  else if (type == "EXCH_M") {v = this->v0.at("EXCH_M"); direct=true; s=this->s12;}
  else {throw psi::PSIEXCEPTION("Incorrect type of matrix element");}

  if (direct) {return this->overlap_corrected_direct  (v   );}
  else        {return this->overlap_corrected_indirect(v, s);}
}
double TIData::overlap_corrected_direct(void) 
{
  return this->overlap_corrected_indirect(0.0, this->s12);
}
double TIData::overlap_corrected_direct(double v) 
{
  return v / (1.0 - this->s12*this->s12);
}
double TIData::overlap_corrected_indirect(double v, double s) 
{
  double E1 = this->e1; 
  double E2 = this->e2;
  if (!this->diagonal_correction) {E1-=this->de1; E2-=this->de2;}
  return 1.0/(1 - s*s) * (v - s*((double)0.5)*(E1 + E2));
}
double TIData::coupling_direct_coul(void) {
  double v;
  // Determine Coulombic coupling
  if (this->trcamm_approximation) {v = this->v0_trcamm->level(this->trcamm_convergence)->get(0,0);}
  else {v = this->v0.at("COUL");}
  return v;
}
double TIData::coupling_direct_exch(void) {
  double v;
  // Determine pure exchange coupling
  if (this->mulliken_approximation) {v = this->v0.at("EXCH_M");}
  else {v = this->v0.at("EXCH");}
  return v;
}
double TIData::coupling_direct(void) {

  double V_Coul = this->coupling_direct_coul();
  double V_Exch = this->coupling_direct_exch();
  double V_Ovrl = 0.0;

  // Do overlap correction?
  if (this->overlap_correction) {
      V_Coul = this->overlap_corrected_direct(V_Coul);
      V_Exch = this->overlap_corrected_direct(V_Exch);
      V_Ovrl = this->overlap_corrected_direct();
  }
  //
  return V_Coul + V_Exch + V_Ovrl;
}

double TIData::coupling_indirect(void) {
  double V_2 = this->coupling_indirect_ti2();
  double V_3 = this->coupling_indirect_ti3();
  return V_2 + V_3;
}
double TIData::coupling_indirect_ti2(void) {
  double E1 = this->e1; 
  double E2 = this->e2;
  if (!this->diagonal_correction) {E1-=this->de1; E2-=this->de2;}
  double E3 = this->e3;
  double E4 = this->e4;
  double V_ET1 = this->v0.at("ET1");
  double V_ET2 = this->v0.at("ET2");
  double V_HT1 = this->v0.at("HT1");
  double V_HT2 = this->v0.at("HT2");

  if (this->overlap_correction) {
      V_ET1= this->overlap_corrected("ET1");
      V_ET2= this->overlap_corrected("ET2");
      V_HT1= this->overlap_corrected("HT1");
      V_HT2= this->overlap_corrected("HT2");
  }

  double v =-(V_ET1*V_HT2)/(E3-E1) -(V_HT1*V_ET2)/(E4-E1);
  return v;
}
double TIData::coupling_indirect_ti3(void) {
  double E1 = this->e1; 
  double E2 = this->e2;
  if (!this->diagonal_correction) {E1-=this->de1; E2-=this->de2;}
  double E3 = this->e3;
  double E4 = this->e4;
  double V_ET1 = this->v0.at("ET1");
  double V_ET2 = this->v0.at("ET2");
  double V_HT1 = this->v0.at("HT1");
  double V_HT2 = this->v0.at("HT2");
  double V_CT;
  if (!this->mulliken_approximation) 
     {V_CT  = this->v0.at("CT");}
  else 
     {V_CT  = this->v0.at("CT_M");}

  if (this->overlap_correction) {
      V_ET1= this->overlap_corrected("ET1");
      V_ET2= this->overlap_corrected("ET2");
      V_HT1= this->overlap_corrected("HT1");
      V_HT2= this->overlap_corrected("HT2");
      V_CT = this->overlap_corrected_indirect(V_CT, this->s34);
  }

  double v = (V_ET1*V_ET2 + V_HT1*V_HT2) * V_CT / ((E3-E1)*(E4-E1));
  return v;
}
double TIData::coupling_total(void) {
 return this->coupling_direct() + this->coupling_indirect();
}
