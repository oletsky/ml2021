import numpy as np;
#Solving linear algebraic equations
matr=[
        [1,1],
        [1, 1]
        
     ];
v=[10,2];
x=np.linalg.solve(matr,v);
print (x);
