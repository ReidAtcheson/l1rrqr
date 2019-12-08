# l1rrqr
L1-norm based rank-revealing factorization


This implementation is based off my paper [A Rank Revealing Factorization Using Arbitrary Norms](https://arxiv.org/abs/1905.02355).


I give the option to use a NAG library solver for the resulting least-1-norm solutions because it is faster and more accurate than
using a generic linear program solver, but if the module is not found it will default to using the linear program solver in scipy
instead.



# Explanation of supplementary files

## Condition number tests
The file `condition_number_tests.py` shows that the `Q` factor in `AP=QR` has `1-norm` condition number independent of the
`1-norm` condition number of `A`. This appears to be largely true, but when the conditioning of `A` becomes very large
the `l1rrqr` algorithm starts to fail. This is likely because it does not have logic for handling the rank deficient case,
which is fixable.


## Factorization error tests

The file `factorization_error_tests.py` shows that the error `AP-QR` measured in the `1-norm` is independent of the 
`1-norm` condition number of `A`. Again the factorization appears to fail in the rank deficient case, so in these
cases the error is very large


## Nonincreasing diagonal

The file `nonincreasing_diagonal.py` prints the diagonal of `R` for a factorization of a random matrix. 
