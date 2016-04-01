#coding:utf-8
#
#  function [f, J] = fg3dcircle(a, X, w) 
#% --------------------------------------------------------------------- 
#% FG3DCIRCLE.M   Function and gradient calculation for  
#%                least-squares circle fit. 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% --------------------------------------------------------------------- 
#% Input  
#% a        Parameters [x0 y0 z0 alpha beta s]'. 
#%          Dimension: 6 x 1. 
#% 
#% X        Array [x y z] where x = vector of x-coordinates,  
#%          y = vector of y-coordinates and z = vector of z-coordinates.  
#%          Dimension: m x 3.  
#% 
#%  
#%  
#% Output 
#% f        Vector [f1; f2] where f1 = vector of distances from the  
#%          points to the plane containing the circle and f2 = vector  
#%          of distances from the points to the cylinder containing  
#%          the circle.  
#%          Dimension: 2m x 1.  
#% 
#%   
#% 
#% Modular structure: FGRROT3.M, FRROT3.M, DRROT3.M.  
#% 
#% [f <, J >] = fg3dcircle(a, X <, w >) 
#% --------------------------------------------------------------------- 

import numpy as np
from fgrrot3 import *
from frrot3 import *
from drrot3 import *


def fg3dcircle(a=None, X=None, w=None): 
    """FG3DCIRCLE.py   Function and gradient calculation for least-squares circle fit. 
    Input  
    a        Parameters [x0 y0 z0 alpha beta s]'. 
              Dimension: 6 x 1. 
              
    X        Array [x y z] where x = vector of x-coordinates,  
              y = vector of y-coordinates and z = vector of z-coordinates.  
              Dimension: m x 3.  
    Output 
    f        Vector [f1; f2] where f1 = vector of distances from the  
              points to the plane containing the circle and f2 = vector  
            of distances from the points to the cylinder containing  the circle.  
            Dimension: 2m x 1.  
     Modular structure: FGRROT3.py, FRROT3.py, DRROT3.py.  
    [f <, J >] = fg3dcircle(a, X <, w >) """
    m = np.size(X, 0);  
    M1=np.mat(np.ones((m, 1)))
    MZ=np.mat(np.zeros((m, 1)))
    #if no weights are specified, use unit weights  
    if w==None:
        w = M1; 
    if np.size(a,1)!=1:
        raise ArithmeticError,"a is not m*1"   
    #µÃµ½x0 y0 z0
    x0 = a[0,0];
    y0 = a[1,0];
    z0 = a[2,0]; 

    alpha = a[3,0];
    beta = a[4,0];

    s = a[5,0];
    
    R, DR1,DR2,DR3 = fgrrot3(theta=np.mat([[alpha], [beta], [0]]));
    #Î»ÒÆ¾ØÕó    
    Xt = (X -M1  * [x0,y0,z0] )* R.T;  
    xt = np.asmatrix(Xt[:, 0]);
    yt = np.asmatrix(Xt[:, 1]);
    zt = np.asmatrix(Xt[:, 2]);  
    
    rt = np.sqrt(np.multiply(xt,xt) + np.multiply(yt,yt));  
    Nt1 = np.concatenate(((xt/rt),(yt/rt),MZ),1); 
    f1 = np.sum(np.multiply(Xt, Nt1),1)#np.dot(Xt, Nt1, 2);
    f1 = f1 - s;  
    #¼ÓÈ¨
    f1 = np.multiply(w,f1); # incorporate weights 
     
    Nt2 = np.concatenate((np.zeros((m, 2)),M1),1);  
    f2 = np.sum(np.multiply(Xt[:,2], Nt2),1)#np.dot(Xt[:,2], Nt2, 2);  
    #¼ÓÈ¨
    f2 = np.multiply(w,f2); # incorporate weights  
    f =np.concatenate((f1, f2),0);
       
    #if nargout > 1:
    #form the Jacobian matrix
    J1 = np.mat(np.zeros((m, 6))); #derivatives of f1  
    J2 = np.mat(np.zeros((m, 6))); #derivatives of f2  
    A1 = M1* (R * np.mat("-1;0;0")).T;  
    #print np.sum(np.multiply(A1, Nt1),1)
    J1[:, 0] = np.sum(np.multiply(A1, Nt1),1)#np.dot(A1, Nt1, 2);  
    J2[:, 0] = np.sum(np.multiply(A1, Nt2),1)#np.dot(A1, Nt2, 2);  
    A2 = M1 * (R * np.mat("0;-1;0")).T;  
    J1[:, 1] = np.sum(np.multiply(A2, Nt1),1)#np.dot(A2, Nt1, 2); 
    J2[:, 1] = np.sum(np.multiply(A2, Nt2),1)#np.dot(A2, Nt2, 2);  
    A3 = M1* (R * np.mat("0 ;0 ;-1")).T; 
    J1[:, 2] = np.sum(np.multiply(A3, Nt1),1)#np.dot(A3, Nt1, 2);  
    J2[:, 2] = np.sum(np.multiply(A3, Nt2),1)#np.dot(A3, Nt2, 2);  
    A4 = (X - M1 * np.mat([x0,y0,z0])) * DR1.T;  
    J1[:, 3] = np.sum(np.multiply(A4, Nt1),1)#np.dot(A4, Nt1, 2);  
    J2[:, 3] = np.sum(np.multiply(A4, Nt2),1)#np.dot(A4, Nt2, 2); 
    A5 = (X -M1 * np.mat([x0,y0,z0])) * DR2.T;  
    J1[:, 4] =np.sum(np.multiply(A5, Nt1),1)# np.dot(A5, Nt1, 2);  
    J2[:, 4] = np.sum(np.multiply(A5, Nt2),1)#np.dot(A5, Nt2, 2); 
    J1[:, 5] =-1*(np.ones((m, 1)));  
    J2[:, 5] = MZ;  
    W = np.mat(np.diag(np.array(w[:,0])[:,0]));  
    J1 = W * J1; # incorporate weights  
    J2 = W * J2; # incorporate weights  
    J = np.concatenate((J1,J2),0); 
    return f, J 
#  end % if nargout 
if __name__=="__main__":
    import math
    print fg3dcircle.__doc__
    a=np.mat("1 2 3;4 5 6;7 8 9")
    b=np.mat("   8     1     6;3     5     7;4     9     2")
    #print np.dot(a,b.T)
    #print np.tensordot(a,b,0)
    theta=np.mat([[ 1.40048515],
                  [ 2.10329182],
                  [ 1.12845222]]);
    c11=np.mat([[182.86779073979     ,2325.87529134155     ,-294.40450739032 ],   
[160.70160950720     ,2406.12086429678    , -294.58945521292] ,
[213.76616244049    , 2504.70831446148    , -294.68201477650],
[421.01919061425    , 2377.73490680325     ,-294.38556654694],
[365.44126911669    , 2290.78641632270     ,-294.41985361918 ],
[271.79572682184    , 2269.46402226983    , -294.34976761283 ],
[222.96860277392    , 2287.43892234611     ,-294.44470758825 ]])
    #print theta
    a=np.mat([291.63283003918434,2399.1548675538656,-294.50802472282004,math.pi/180*5,math.pi/180*80,0])
    f,j = fg3dcircle(a.T,c11 )
    print "f=",f,"\nj=",j



