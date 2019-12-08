import numpy as np
import scipy.linalg as la
from l1rrqr import l1rrqr



#These tests show that the factorization error AP-QR for the l1rrqr algorithm
#is independent of the conditioning of A

seed=234234
np.random.seed(seed)

nconds=10
minorder=5
maxorder=30
print("Matrix order,Condition of full matrix,||AP-QR||_1")
for m in range(minorder,maxorder):
    for cond in np.logspace(0,nconds,m):
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
        print("{},{},{}".format(m,np.linalg.cond(M,p=1),np.linalg.norm(M[:,p]-Q@R,ord=1)))
