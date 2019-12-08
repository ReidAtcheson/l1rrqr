# l1rrqr
L1-norm based rank-revealing factorization


This implementation is based off my paper [A Rank Revealing Factorization Using Arbitrary Norms](https://arxiv.org/abs/1905.02355).


I give the option to use a NAG library solver for the resulting least-1-norm solutions because it is faster and more accurate than
using a generic linear program solver, but if the module is not found it will default to using the linear program solver in scipy
instead.
