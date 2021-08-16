import random as rand

def counter(mas):
    count=dict()
    n=len(mas)
    for i in range(n):
        v=mas[i]
        if v in count:
            count[v]=count[v]+1
        else:
            count[v]=1
    return count


def getrandbyprobs(p):
    eps=1.e-10
    sm=sum(p)
    if (sm<1.-eps) or (sm>1.+eps):
        raise ValueError("Sum ",sm," not equals 1!")
    n=len(p)
    r=rand.random()
    s=0
    ind=-1
    for i in range(n):
        ind=i
        if p[ind]<eps:
            continue
        s=s+p[ind]
        if r<s:
            break
    return ind

p=[0.6,0.3,0.1]

#Checking
k=100000
mas=[]
for i in range(k):
    rr=getrandbyprobs(p)
    mas.append(rr)
print("Experimental results:")
ll=len(mas)
counts=counter(mas)
for k in counts.keys():
    print(k,"-",counts[k]," - ",counts[k]/ll)

    
    
