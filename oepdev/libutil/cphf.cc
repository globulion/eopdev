#include "cphf.h"

namespace oepdev{

using namespace std;
using namespace psi;


CPHF::CPHF(SharedWavefunction ref_wfn, Options& options) :
            _primary(ref_wfn->basisset()),
            _cocc(ref_wfn->Ca_subset("AO","OCC")),
            _cvir(ref_wfn->Ca_subset("AO","VIR")),
            _eps_occ(ref_wfn->epsilon_a_subset("AO","OCC")),
            _eps_vir(ref_wfn->epsilon_a_subset("AO","VIR")),
            _memory(Process::environment.get_memory()),
            _options(options),
            _diis_dim(options.get_int("CPHF_DIIS_DIM")),
            _no(ref_wfn->doccpi()[0]),
            _nv(ref_wfn->nmo()-ref_wfn->doccpi()[0]),
            _nn(ref_wfn->basisset()->nbf())
{
    if (not (_options.get_str("REFERENCE") == (std::string)"RHF")) {
       throw PSIEXCEPTION("oepdev_util::CPHF now only for RHF wavefunctions!");
    }
    _maxiter   = _options.get_int("CPHF_MAXITER");
    _conv      = _options.get_double("CPHF_CONVERGENCE");
    _with_diis = _options.get_bool("CPHF_DIIS");

    // Initialize DIIS manager
    if (_with_diis) {
        for (unsigned int z=0; z<3; z++) {
             #if OEPDEV_USE_PSI4_DIIS_MANAGER == 0
               _diis.push_back(std::shared_ptr<oepdev::DIISManager>(new oepdev::DIISManager(_diis_dim, _no, _nv)));
             #else
               _diis.push_back(std::shared_ptr<psi::DIISManager>(new psi::DIISManager(_diis_dim, "CPHF DIIS", 
                                 psi::DIISManager::LargestError, psi::DIISManager::InCore)));
               _diis[z]->set_error_vector_size(1, psi::DIISEntry::Matrix, 
                          std::shared_ptr<Matrix>(new Matrix("Error Vector", _no, _nv)).get());
               _diis[z]->set_vector_size(1, psi::DIISEntry::Matrix, 
                          std::shared_ptr<Matrix>(new Matrix("DIIS Vector", _no, _nv)).get());
             #endif
        }
    }
}

CPHF::~CPHF() 
{
}

void CPHF::compute(void) {

    // Setup dipole perturbations //

    std::vector<std::shared_ptr<Matrix>> Fmo;
    std::vector<std::shared_ptr<Matrix>> Fao;

    for (unsigned int z = 0; z < 3; z++) {
      Fao.push_back(std::shared_ptr<Matrix>(new Matrix("F", _nn, _nn)));
    }
   
    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(_primary)); 
    std::shared_ptr<OneBodyAOInt> ints(fact->ao_dipole());

    ints->compute(Fao);

    for (unsigned int z=0; z<3; z++) {
      Fmo.push_back(Matrix::triplet(_cocc,Fao[z],_cvir,true,false,false));
      Fmo[z]->scale(-2.0);
    }

    // Setup uncoupled guess //

    double* eop = _eps_occ->pointer();
    double* evp = _eps_vir->pointer();

    std::vector<std::shared_ptr<Matrix>> Xmo;
    
    for (unsigned int z=0; z<3; z++) {
       Xmo.push_back(std::shared_ptr<Matrix>(new Matrix("X",_no,_nv)));
       Xmo[z]->copy(Fmo[z]);
       double** Xp = Xmo[z]->pointer();
       for (unsigned int i=0; i<_no; i++) {
            for (unsigned int a=0; a<_nv; a++) {
                 Xp[i][a] /= (evp[a] - eop[i]);
            }
       }
    }

    // Setup history for convergence //

    std::vector<std::shared_ptr<Matrix>> Xmo_old;
    for (unsigned int z=0; z<3; z++) {
       Xmo_old.push_back(std::shared_ptr<Matrix>(new Matrix("X",_no,_nv)));
    }
 
    // Set JK object for 2-body part //

    std::shared_ptr<JK> jk = JK::build_JK(_primary, BasisSet::zero_ao_basis_set(), _options);
    jk->set_memory(_memory / 8L); 
    jk->initialize();
    jk->print_header();
    
    std::vector<std::shared_ptr<Matrix>>& Cl = jk->C_left();
    std::vector<std::shared_ptr<Matrix>>& Cr = jk->C_right();

    const std::vector<std::shared_ptr<Matrix>>& J = jk->J();
    const std::vector<std::shared_ptr<Matrix>>& K = jk->K();

    bool converged = false;
    double max_rms = 0.0;
    outfile->Printf("  => CPHF Iterations <=\n\n");
  
    for (unsigned int iter=0; iter < _maxiter; iter++) {

        // Compute generalized Fock matrices

        Cl.clear();
        Cr.clear();
        for (unsigned int z=0; z<3; z++) {
             Cl.push_back(_cocc);
             Cr.push_back(Matrix::doublet(_cvir,Xmo[z],false,true));
        }

        jk->compute();

        // Update amplitudes

        for (unsigned int z=0; z<3; z++) {
             Xmo[z]->copy(Fmo[z]);
             std::shared_ptr<Matrix> J1 = Matrix::triplet(
                 _cocc,J[z],_cvir,true,false,false);
             std::shared_ptr<Matrix> K1 = Matrix::triplet(
                 _cocc,K[z],_cvir,true,true ,false);
             std::shared_ptr<Matrix> K2 = Matrix::triplet(
                 _cocc,K[z],_cvir,true,false,false);
             J1->scale( 4.0);
             K1->scale(-1.0);
             K2->scale(-1.0);
             Xmo[z]->subtract(J1);
             Xmo[z]->subtract(K1);
             Xmo[z]->subtract(K2);
             double** Xp = Xmo[z]->pointer();
             for (unsigned int i=0; i<_no; i++) {
                 for (unsigned int a=0; a<_nv; a++) {
                     Xp[i][a] /= (evp[a] - eop[i]);
                 }
             }

             // use DIIS
             Xmo_old[z]->subtract(Xmo[z]);
             if (_with_diis) {
                 #if OEPDEV_USE_PSI4_DIIS_MANAGER == 0
                   _diis[z]->put(Xmo_old[z], Xmo[z]);  
                   _diis[z]->compute();
                   _diis[z]->update(Xmo[z]); 
                 #else
                   _diis[z]->add_entry(2, Xmo_old[z].get(), Xmo[z].get());
                   _diis[z]->extrapolate(1, Xmo[z].get());
                 #endif
              }
        }

        // Compute amplitude convergence

        max_rms = 0.0;
        for (unsigned int z=0; z<3; z++) {
             double rms = Xmo_old[z]->rms();
             max_rms = std::max(rms, max_rms);
             Xmo_old[z]->copy(Xmo[z]);
        }

        outfile->Printf("  @CPHF Iter %3d: %11.3E\n", iter + 1, max_rms);
     
        if (max_rms < _conv) {  
            converged = true; 
            break;
        }
    } 

    if (converged) {
        outfile->Printf("\n CPHF Converged.\n\n");
    } else {
        outfile->Printf("\n CPHF Failed.\n\n");
    }

    // Compute and print the dipole polarizability tensor
    
    std::shared_ptr<Matrix> P(new Matrix("Polarizability", 3, 3));
    double** Pp = P->pointer();
    for (unsigned int z1=0; z1<3; z1++) {
        for (unsigned int z2=0; z2<3; z2++) {
            Pp[z1][z2] = Xmo[z1]->vector_dot(Fmo[z2]);
        }
    }
    P->set_name("Molecular Polarizability");
    _molecular_polarizability = P;

}

void CPHF::print(void) const {
    outfile->Printf("  => Setup <=\n\n");
    outfile->Printf("   Memory      = %11zu MB\n" , _memory / (1024L*1024L));
    outfile->Printf("   Maxiter     = %11d\n"     , _maxiter);
    outfile->Printf("   Convergence = %11.3E\n"   , _conv);
    outfile->Printf("   Nocc        = %11d\n"     , _no);
    outfile->Printf("   Nvir        = %11d\n"     , _nv);
    outfile->Printf("   Nbasis      = %11d\n"     , _nn);
    outfile->Printf("\n");
    outfile->Printf("  => Molecule <=\n\n");
   _primary->molecule()->print();
    outfile->Printf("  => Primary Basis <=\n\n");
   _primary->print();
}

}// EndNameSpace oepdev
