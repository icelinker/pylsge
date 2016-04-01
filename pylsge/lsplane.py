#tested beta
#  function [x0, a, d, normd] = lsplane(X) 
#% --------------------------------------------------------------------- 
#% LSPLANE.M   Least-squares plane (orthogonal distance 
#%             regression). 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% --------------------------------------------------------------------- 
#% Input     
#% X        Array [x y z] where x = vector of x-coordinates,  
#%          y = vector of y-coordinates and z = vector of  
#%          z-coordinates.  
#%          Dimension: m x 3.  
#%  
#% Output    
#% x0       Centroid of the data = point on the best-fit plane. 
#%          Dimension: 3 x 1.  
#%  
#% a        Direction cosines of the normal to the best-fit  
#%          plane.  
#%          Dimension: 3 x 1. 
#%  
#%  
#% 
#% [x0, a <, d, normd >] = lsplane(X) 
#% --------------------------------------------------------------------- 


import numpy as np
from scipy import linalg
def lsplane(X):
    """function [x0, a, d, normd] = lsplane(X)
    Input
    X       Array [x y z] where x = vector of x-coordinates,  
            y = vector of y-coordinates and z = vector of  
            z-coordinates.  
            Dimension: m x 3.  
  
    Output    
    x0      Centroid of the data = point on the best-fit plane. 
            Dimension: 3 x 1.  
  
    a       Direction cosines of the normal to the best-fit
            plane.  
            Dimension: 3 x 1. 
    [x0, a <, d, normd >] = lsplane(X) """
    #check number of data points 
    if isinstance(X,np.matrix):
        x=X;
    else:
        x=np.mat(X);
        
    m = x.shape[0]; 
    if m < 3:
        raise ArithmeticError,"At least 3 data points required: "
    else:
        #calculate centroid 
        x0 = np.mean(x,0).T; 
    #form matrix A of translated points 
    A = np.concatenate(((x[:,0] - x0[0]),(x[:,1] - x0[1]),(x[:,2] - x0[2])),1);
    #calculate the SVD of A 
    #find the smallest singular value in S and extract from V the 
    #corresponding right singular vector
    #[U,S,V] = svd(X,0) produces the "economy size" decomposition. If X is m-by-n with m > n, then svd computes only the first n columns of U and S is n-by-n. 
    try:
        U, S, V =np.linalg.svd(A,0)
        #print 'U:',U,"\nS=",S,'\nV=',V
    except :
        raise ArithmeticError, "SVD fail"
        
        
    #find the smallest singular value in S and extract from V the 
    #corresponding right singular vector
    s= S.min(); 
    #print 's=',s
    i=S.argmin()
    #print "i=",i
    a = V.T[:, i] 
    #calculate residual distances, if required   
    d = U[:, i]*s; 
    normd = np.linalg.norm(d); 
    return x0, a, d, normd
if __name__=="__main__":
    print lsplane.__doc__
    #print lsplane.__name__
    c11=np.mat([[277.804,    -2736.349,  -452.272 ],   
         [867.145,    -2660.725,  -522.91] ,
         [1031.517,   -2694.394,  -633.68],
         [977.183,    -2680.959,  -593.538]])
    #print c11
    x0, a, d, normd=lsplane(c11)
    #print 'x0=',x0,'\na=',a,'\nd=',d,'\nnormd=',normd
    





