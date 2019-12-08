# l1rrqr
L1-norm based rank-revealing factorization


This implementation is based off my paper [A Rank Revealing Factorization Using Arbitrary Norms](https://arxiv.org/abs/1905.02355).


I give the option to use a NAG library solver for the resulting least-1-norm solutions because it is faster and more accurate than
using a generic linear program solver, but if the module is not found it will default to using the linear program solver in scipy
instead.



# Explanation of supplementary files

## Condition number tests
The file `condition_number_tests.py` shows that the `Q` factor in `AP=QR` has `1-norm` condition number independent of the
`1-norm` condition number of `A`. This appears to be largely true, but larger condition numbers don't work because `l1rrqr`
returns a reduced factorization for the rank deficient case and so `Q` is rectangular in those cases.

## Factorization error tests

The file `factorization_error_tests.py` shows that the error `AP-QR` measured in the `1-norm` is independent of the 
`1-norm` condition number of `A`. 

## Nonincreasing diagonal

The file `nonincreasing_diagonal.py` prints the diagonal of `R` for a factorization of a random matrix. 
