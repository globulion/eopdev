void scf_perturb(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt)
{
   // Set up HF superfunctional
   std::shared_ptr<psi::SuperFunctional> func = oepdev::create_superfunctional("HF", opt);

   // Initialize the perturbed wavefunction
   std::shared_ptr<oepdev::RHFPerturbed> scf = std::make_shared<oepdev::RHFPerturbed>(wfn, func, opt, wfn->psio());

   /* Perturb the system with the uniform electric field [Fx, Fy, Fz].
      Then, add two point charges of charge qi placed at [Rxi, Ryi, Rzi].
      Provide all these values in atomic units! */
   const double Fx = 0.04, Fy = 0.05, Fz = -0.09;
   const double Rx1= 0.00, Rx2= 1.30, Rx3= -1.00;
   const double Rx1= 0.10, Rx2=-0.30, Rx3=  3.50;
   const double q1 = 0.30, q2 =-0.09; 

   scf->set_perturbation(Fx, Fy, Fz);        /* set it only once, setting it again will overwrite the field, not add */
   scf->set_perturbation(Rx1, Ry1, Rz1, q1); 
   scf->set_perturbation(Rx2, Ry2, Rz2, q2); /* more charges can be added */

   // Solve perturbed SCF equations 
   scf->compute_energy();

   // Grab some data
   double energy = scf->reference_energy();      // Total energy of the system
   std::shared_ptr<psi::Matrix> Da = scf->Da();  // One-particle density matrix

   /* Note that the external field and charges perturb only one-electron Hamiltonian.*/
}
