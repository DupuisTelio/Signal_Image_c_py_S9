import numpy as np

def root(t,T):
    r = t
    while T[r] !=r:
        r = T[r]
    return r

def label(I):
    dim_y,dim_x=np.shape(I)
    Ia = np.zeros((dim_y+2,dim_x+2),dtype=int)
    Ia[1:1+dim_y,1:1+dim_x]=I
    T = [0]
    cpt = 1
    for y in range (1,1+dim_y):
        for x in range (1,1+dim_x):
            if Ia[y,x]:
                t=[Ia[y-1,x-1],Ia[y-1,x],Ia[y-1,x+1],Ia[y,x-1]]
                m=0
                for k in range(0,4):
                    if t[k]:
                        t[k] = root(t[k],T)
                        if m ==0:
                            m=t[k]
                        else:
                            m=min(m,t[k])
                if m==0:
                    Ia[y,x] = cpt
                    T.append(cpt)
                    cpt=cpt+1
                else:
                    Ia[y,x] = m
                    for k in range(0,4):
                        if t[k] and t[k] !=m:
                            T[t[k]]=m

    print(cpt)
    cpt=0
    for i in range(0,len(T)):
        if T[i] == i:
            T[i] = cpt
            cpt= cpt+1
        else:
            T[i]=T[T[i]]
    print(cpt)
    for y in range(1,1+dim_y):
        for x in range(1,1+dim_x):
            Ia[y,x]=T[Ia[y,x]]
    return Ia[1:1+dim_y, 1:1+dim_x]

