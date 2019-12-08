import numpy as np
import scipy.linalg as la
from l1rrqr import l1rrqr


#These tests demonstrate that while the `Q` in `AP=QR` from the `l1rrqr` algorithm
#is not orthgonal the way it is in a traditional rank-revealing QR, its condition
#number is still independent of the conditioning of A.

seed=234234
np.random.seed(seed)

nconds=7
minorder=5
maxorder=30
print("Matrix order,Condition of full matrix,Condition of Q")
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
        print("{},{},{}".format(m,np.linalg.cond(M,p=1),np.linalg.cond(Q,p=1)))
