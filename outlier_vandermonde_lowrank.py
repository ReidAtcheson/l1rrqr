import numpy as np
import scipy.linalg as la
from l1rrqr import l1rrqr


m=15
n=15
k=3
V=np.zeros((m,n))
noise=np.zeros((m,n))
noise[5,5]=1
xs=np.linspace(-1.0,1.0,m)
V[:,0]=1.0
for i in range(1,n):
    V[:,i]=V[:,i-1]*xs



Q1,R1,p1=la.qr(V+noise,pivoting = True)
P1=np.identity(m)[:,p1]
Q2,R2,p2=l1rrqr(V+noise)
P2=np.identity(m)[:,p2]


print("               L2-RRQR,                        L1-RRQR")
print("L1-error:      {},                {}".format(np.linalg.norm(V[:,p1]-Q1[:,0:k]@R1[0:k,:],ord=1),np.linalg.norm(V[:,p2]-Q2[:,0:k]@R2[0:k,:],ord=1)))
print("L2-error:      {},                {}".format(np.linalg.norm(V[:,p1]-Q1[:,0:k]@R1[0:k,:],ord=2),np.linalg.norm(V[:,p2]-Q2[:,0:k]@R2[0:k,:],ord=2)))
