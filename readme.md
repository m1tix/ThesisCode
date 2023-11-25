# Code for Thesis

Public repository for code related to my thesis "The Hasse Norm Theorem and Norm
Equations". Everything is written in Python and requires Sage. One can run a
specific file in the Sage interpreter with 'load()' and 'attach()'

The following files are included:

- [biquadratic_norm](biquadratic_norm.sage): code to generate Examples 7.2.8 and
  7.2.9.
- [counterexample_Sunits](counterexample_Sunits.sage): code to generate the
  counterexamples for the suitability of certain subsets of non-Galois extensions
  (remark to Lemma 7.2.7).
- [fincke_pohst](fincke_pohst.sage): code implementing the Fincke-Pohst
  algorithm to solve absolute integral norm equations.
- [local_everywhere_quadratic](local_everywhere_quadratic.sage): code to generate a quadratic extension
  and a prime number q such that there exists no integral element of norm q, but
  the associated norm equation is locally soluble everywhere.
- [norm_example_baker](norm_example_baker.sage): sloppy code implementing
  Baker's method for solving integral norm equations.
- [ploty](ploty.py): code for the plots of Figure 1 using data from
  [here](data).
