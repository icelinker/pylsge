#  function [x0, a, d, normd] = ls3dline(X)
#% ---------------------------------------------------------------------
#% LS3DLINE.M   Least-squares line in 3 dimensions.
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
#% x0       Centroid of the data = point on the best-fit line.
#%          Dimension: 3 x 1. 
#% 
#% a        Direction cosines of the best-fit line. 
#%          Dimension: 3 x 1.
#% 
#% <Optional... 
#% d        Residuals. 
#%          Dimension: m x 1. 
#% 
#% normd    Norm of residual errors. 
#%          Dimension: 1 x 1. 
#% ...>
#%
#% [x0, a <, d, normd >] = ls3dline(X)
#% ---------------------------------------------------------------------
import numpy as np
def ls3dline(X=None):
    """Input
    X   Array [x y z] where x = vector of x-coordinates,
        y = vector of y-coordinates and z = vector of 
        z-coordinates. 
        Dimension: m x 3.
    Output
    x0  Centroid of the data = point on the best-fit line.
        Dimension: 3 x 1.  
    a   Direction cosines of the best-fit line. 
        Dimension: 3 x 1. 
    d   Residuals. 
        Dimension: m x 1. 
    normd    Norm of residual errors. 
            Dimension: 1 x 1. 
    [x0, a <, d, normd >] = ls3dline(X)"""
    #check number of data points 
    if X==None:
        raise ArithmeticError,'Require [x y z] array'
    else:
        if isinstance(X,np.matrix):
            pass
        else:
            x=np.mat(X);  
        m = np.size(X, 0);
        if m < 3:
            raise ArithmeticError,'At least 3 data points required: ' 
    #end
    #calculate centroid
    x0 = np.mat(np.mean(X,0).T);
    #form matrix A of translated points
    A = np.concatenate(((X[:, 0] - x0[0,0]),(X[:, 1] - x0[1,0]),(X[:, 2] - x0[2,0])),1);
    #calculate the SVD of A
    U, S, V = np.linalg.svd(A, 0);
    #find the largest singular value in S and extract from V the
    #corresponding right singular vector
    s = S.max();
    i=S.argmax()
    a = V.T[:, i];
    #calculate residual distances, if required  
    m = np.size(X, 0); 
    d = np.mat(np.zeros((m, 1))); 
    for i in range(0,m):
        dd=X[i, 0:3]- x0.T
        dc=np.cross(dd, a.T)
        d[i] = np.linalg.norm(dc); 
    #end % for i 
    
    normd = np.linalg.norm(d); 
    #end % if nargout 
    return  x0, a , d, normd 
if __name__=="__main__":
    print ls3dline.__doc__
    x=np.mat("1.001,1.0,1.002;2.0,2.005,2.009;3.008,3.006,3.005;4.001,4.0,4.0")
    x0, a , d, normd =ls3dline(x)
    print "x0=", type(x0),x0,'\na=', type(a),a , '\nd=', type(d), d, '\nnormd=', type(normd),normd
#x0 =
#
#    2.5025
#    2.5027
#    2.5040
#
#
#a =
#
#    0.5778
#    0.5774
#    0.5768
#
#
#d =
#
#    0.0023
#    0.0048
#    0.0025
#    0.0003
#
#
#normd =
#
#    0.0059

