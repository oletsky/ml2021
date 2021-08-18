import numpy as np;
#Solving linear algebraic equations by psevdoinversion
matr=[
        [1],
        [3]
        
    ]
v=[1, 9]

pv=np.linalg.pinv(matr);
print("Pseudoinverse matrix:")
print (pv)
x=np.dot(pv,v)
print("Resulting vector:")
print(x)
