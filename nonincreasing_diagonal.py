import numpy as np
import scipy.linalg as la
from l1rrqr import l1rrqr



#This demonstrates the diagonal of `R` is nonincreasing,
#the key property for a rank revealing factorization

seed=234234
np.random.seed(seed)

maxorder=30
m=30
cond=1e7

#Generate a random matrix with 
#prescribed conditioning by working
#backwards from an SVD
A=np.random.rand(m,m)
B=np.random.rand(m,m)
U,_=np.linalg.qr(A)
V,_=np.linalg.qr(B)
D=np.diag(np.linspace(1.0,cond,m))
M=U@D@V
Q,R,p=l1rrqr(M)
print(np.diag(R))
