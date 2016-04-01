import numpy as np
def crosspoint(v,vp,n,np):
    if v.shape==(3,1)and n.shape==(3,1)and vp.shape==(3,1)and np.shape==(3,1):
        pass
    else:
        return (np.mat("0;0;0"),False)
    dt=-n.T*(vp-np)/(v.T*n)
    dp=vp+v*dt
    print dp
    return dp
if __name__=="__main__":
    v=np.mat("0;0;1")
    vp=np.mat("0;0;1")
    n=np.mat("0;1;1")
    np=np.mat("100;10;-50")
    print crosspoint(v,vp,n,np)
 
