#ifndef _oepdev_libutil_esp_h
#define _oepdev_libutil_esp_h

#include<string>

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "../liboep/oep.h"


namespace oepdev{

using namespace psi;
using namespace std;


using SharedWavefunction       = std::shared_ptr<Wavefunction>;
using SharedMatrix             = std::shared_ptr<Matrix>;
using SharedOEPotential        = std::shared_ptr<OEPotential>;


/** \brief 3-dimensional scalar potential manifold: Abstract base.
 *
 */
class Potential3D
{
  public:
    Potential3D(SharedMatrix pot);
    Potential3D(SharedWavefunction wfn);
    Potential3D(SharedOEPotential oep, const std::string& oepType);

    static std::shared_ptr<Potential3D> build(SharedMatrix pot);
    static std::shared_ptr<Potential3D> build(const std::string& potType, SharedWavefunction wfn);
    static std::shared_ptr<Potential3D> build(SharedOEPotential oep, const std::string& oepType);

    virtual ~Potential3D();

    virtual SharedMatrix get(void) const {return pot_;}

    virtual void compute() = 0;

  protected:
    SharedMatrix pot_;
  private:
    void common_init();
};

class BarePotential3D : public Potential3D
{
  public:
    BarePotential3D(SharedMatrix pot);
    virtual void compute();
  protected:
  private:
    void common_init();
};

class ElectrostaticPotential3D : public Potential3D
{
  public:
    ElectrostaticPotential3D(SharedWavefunction wfn);
    virtual void compute();
  protected:
  private:
    void common_init();
};

class OEPotential3D : public Potential3D
{
  public:
    OEPotential3D(SharedOEPotential oep, const std::string& oepType);
    virtual void compute();
  protected:
  private:
    void common_init();
};

/** \brief Charges from Electrostatic Potential (ESP): abstract base of ESP method.
 *
 */
//class ESP
//{
//  public:
//    //@{ Use atoms as ESP centres
//    ESP(const Wavefunction& wfn);
//    ESP(const OEPotential& oep, const std::string& oepType);
//
//    static std::shared_ptr<ESP> build(const std::string& method, const Wavefunction& wfn);
//    static std::shared_ptr<ESP> build(const std::string& method, const OEPotential& oep, const std::string& oepType);
//    //@}
//
//    virtual ~ESP();
//
//    virtual void fit();
//
//    virtual void compute_potential() = 0;
//    virtual void print_header() const = 0;
//
//  protected:
//    SharedPotential3D potRef_;
//    SharedPotential3D pot_;
//    
//  private:
//    void common_init(void);
//};

//class OEPotentialESP : public ESP
//{
//  public:
//    OEPotentialESP(SharedOEPotential oep, const std::string& oepType);
//  protected:
//  private:
//};

} // EndNameSpace oepdev
#endif //_oepdev_libutil_esp_h
