void iterate(std::shared_ptr<oepdev::IntegralFactory> ints)
{
  // Prepare for direct calculation of ERI's (shell by shell)
  std::shared_ptr<psi::TwoBodyAOInt> tei(ints->eri());

  // Grab the buffer where the integrals for a current shell will be placed
  const double* buffer = tei->buffer();

  // Create iterator to go through all shell quartet combinations
  oepdev::SharedShellsIterator shellIter = oepdev::ShellCombinationsIterator::build(ints, "ALL", 4);

  // Iterate over shells, and then over all integrals in each shell quartet
  for (shellIter->first(); shellIter->is_done() == false; shellIter->next())
  {
       // Compute all integrals between shells in the current quartet
       shellIter->compute_shell(tei);

       // Create iterator to go through all integrals within a shell quartet
       oepdev::SharedAOIntsIterator intsIter = shellIter->ao_iterator("ALL");

       for (intsIter->first(); intsIter->is_done() == false; intsIter->next())
       {
            // Grab current (ij|kl) indices here
            int i = intsIter->i();
            int j = intsIter->j();
            int k = intsIter->k();
            int l = intsIter->l();

            // Grab the (ij|kl) integral
            double integral = buffer[intsIter->index()];
       }
  }
}
