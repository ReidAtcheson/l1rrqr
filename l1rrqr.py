import numpy as np
import scipy.optimize as opt

#Find least-1-norm solution using scipy.optimize
def lst1norm(A,b):
    (m,n)=A.shape
    nvars=m+n
    ncons=2*m
    cons=np.zeros((ncons,nvars))
    cons[0:m,0:m]=-np.identity(m)
    cons[0:m,m:m+n]=A
    cons[m:2*m,0:m]=-np.identity(m)
    cons[m:2*m,m:m+n]=-A
    c=np.zeros(nvars)
    c[0:m]=1.0
    ub=np.zeros(ncons)
    ub[0:m]=b
    ub[m:2*m]=-b
    bounds=[]
    for i in range(0,m):
        bounds.append((0,None))
    for i in range(m,m+n):
        bounds.append((None,None))
    out=opt.linprog(c,cons,ub,None,None,bounds,options={'tol':1e-10,'lstsq' : True})
    return (out.x[m:m+n],out.fun)


#Check to see if NAG module is available,
#if it is then make the default l1-solver 
#to be lst1norm_nag_glin, otherwise
#stick with scipy

#The NAG L1-solver is more robust,accurate, and 
#tends to be a lot faster than the scipy based solver.
l1alg=lst1norm
try:
    from naginterfaces.library.fit import glin_l1sol
    def lst1norm_nag_glin(A,b):
        (m,n)=A.shape
        B=np.zeros((m+2,n+2))
        B[0:m,0:n]=A
        (_,_,soln,objf,rank,i)=glin_l1sol(B,b)
        return (soln[0:n],objf)
    l1alg=lst1norm_nag_glin
except:
    pass



def l1rrqr(A,tol=1e-15):
    (m,n) = A.shape
    Q=np.zeros((m,n))
    R=np.zeros((n,m))
    norms=[np.linalg.norm(A[:,i],ord=1) for i in range(0,n)]
    k=np.argmax(norms)
    perm=[k]
    sout=set(perm)
    sin=set([i for i in range(0,n)]).difference(sout)
    i=0
    Q[:,i]=A[:,k]/np.linalg.norm(A[:,k],ord=1)
    R[i,i]=np.linalg.norm(A[:,k],ord=1)
    while sin:
        i=i+1
        V=Q[:,0:i]
        vals=[l1alg(V,A[:,i])[1] for i in sin]
        ids=list(sin)
        k=ids[np.argmax(vals)]
        sout.add(k)
        sin=sin.difference(sout)
        perm.append(k)
        (c,objf)=l1alg(V,A[:,k])
        print(c)
        print(c.shape)
        print(V.shape)
        y=A[:,k]-V@c
        Q[:,i]=y/np.linalg.norm(y,ord=1)
        R[0:i,i]=c
        R[i,i]=np.linalg.norm(y,ord=1)
    return (Q,R,perm)





