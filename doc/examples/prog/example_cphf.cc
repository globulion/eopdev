void example_cphf(std::shared_ptr<psi::Wavefunction> wfn, psi::Options& opt){

  // build the solver
  std::shared_ptr<oepdev::CPHF> solver = std::make_shared<oepdev::CPHF>(wfn, opt);

  // run the solver to converge CPHF equations
  solver->compute();

  // print the LMO-distributed polarizabilities
  for (int i=0; i<solver->nocc(); i++) {
       solver->polarizability(i)->print();
  }

  // print the molecular polarizability
  solver->polarizability()->print();

  // grab 4th LMO-distributed polarizability and its associated LMO centroid
  psi::SharedMatrix pol_4 = solver->polarizability(3);
  psi::SharedVector rmo_4 = solver->lmo_centroid(3);
};
