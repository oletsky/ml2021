import numpy as np;
#Solving linear algebraic equations by psevdoinversion
matr=[
        [1, 0],
        [0, 1],
        [1, 1]
     ]
m=len(matr)
print("Number of equations: ",m)
n=len(matr[0])
print("Number of variables: ",n)
v=[1, 2, 10]
w=[1./3., 1./3., 1./3.]

#Weighted matrix
wdiag = [ [0]*m for i in range(m) ]
for i in range(m):
    wdiag[i][i]=w[i]

print("Weighting matrix")
print(wdiag)
wmatr=np.dot(wdiag,matr)
wvect=np.dot(wdiag,v)
pv=np.linalg.pinv(wmatr);
print("Pseudoinverse matrix:")
print (pv)
x=np.dot(pv,wvect)
print("Resulting vector:")
print(x)
