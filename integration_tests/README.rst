HelixMC Integration Tests
=========================

This folder contains tests for HelixMC. Currently we only test the helixmc-run
command lines. These test command lines are located in the `cmdline_tests`
folder. To start the test, simply run::

  python runtest.py

Then compare the results with the stored results::

  diff -r test_results stored_results

Code changes that modify the test results should be examined carefully before
commiting. If the changes are expected, you may replace the `stored_results`
folder with the latest `test_results` then commit.

For debugging, one can also just run a particular test. For example::

  python test.py zero_force

This just runs the `zero_force` test in the `cmdline_tests` folder.

To add more tests, just create a runnable command line and put it in a new
file in the `cmdline_tests` folder. Remember to add the `-const_seed` flag to
force deterministic random generator. Do not forget to update the
`stored_results` as well.
