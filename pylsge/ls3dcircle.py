#coding:utf-8
#  function [x0n, an, rn, d, e, f, sigmah, conv, Vx0n, Van, urn, ...  
#            GNlog, a, R0, R] = ls3dcircle(X, x0, a0, r0, tolp, tolg, w) 
#% --------------------------------------------------------------------- 
#% LS3DCIRCLE.M   Least-squares circle in three dimensions using 
#%                Gauss-Newton. 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% --------------------------------------------------------------------- 
#% Input     
#% X        Array [x y z] where x = vector of x-coordinates,  
#%          y = vector of y-coordinates and z = vector of z-coordinates.  
#%          Dimension: m x 3.  
#%  
#% x0       Estimate of the circle centre.  
#%          Dimension: 3 x 1.  
#% 
#% a0       Estimate of the normal to the plane containing the circle.  
#%          Dimension: 3 x 1.  
#%  
#% r0       Estimate of the circle radius.  
#%          Dimension: 1 x 1.  
#%  
#% tolp     Tolerance for test on step length.  
#%          Dimension: 1 x 1.  
#% 
#% tolg     Tolerance for test on gradient. 
#%          Dimension: 1 x 1.  
#%  
#%   
#%  
#% Output   
#% x0n      Estimate of the circle centre.  
#%          Dimension: 3 x 1.  
#%  
#% an       Estimate of the normal direction  
#%          Dimension: 3 x 1.  
#%  
#% rn       Estimate of the circle radius  
#%          Dimension: 1 x 1. 
#%  
#% d        Vector of distances from the points to the circle  
#%          Dimension: m x 1.  
#%  
#% e        Vector of distances from the points to the plane  
#%          containing the circle 
#%          Dimension: m x 1.  
#%  
#% f        Vector of distances from the points to the cylinder  
#%          containing the circle 
#%          Dimension: m x 1.  
#%  
#% sigmah   Estimate of the standard deviation of the weighted  
#%          residual errors.  
#%          Dimension: 1 x 1.  
#%  
#% conv     If conv = 1 the algorithm has converged,  
#%          if conv = 0 the algorithm has not converged 
#%          and x0n, rn, d, e, f and sigmah are current estimates.  
#%          Dimension: 1 x 1.  
#%  
#% Vx0n     Covariance matrix of circle centre.  
#%          Dimension: 3 x 3.  
#% 
#% Van      Covariance matrix of normal direction.  
#%          Dimension: 3 x 3.  
#% 
#% urn      Uncertainty in circle radius.  
#%          Dimension: 1 x 1.  
#%  
#% GNlog    Log of the Gauss-Newton iterations.  
#%          Rows 1 to niter contain  
#%          [iter, norm(f_iter), |step_iter|, |gradient_iter|].  
#%          Row (niter + 1) contains  
#%          [conv, norm(d), 0, 0].  
#%          Dimension: (niter + 1) x 4.  
#%  
#% a        Optimisation parameters at the solution. 
#%          Dimension: 6 x 1.  
#%  
#% R0       Fixed rotation matrix.  
#%          Dimension: 3 x 3.  
#%  
#% R        Upper-triangular factor of the Jacobian matrix 
#%          at the solution.  
#%          Dimension: 6 x 6.  
#% 
#% Modular structure: NLSS11.M, GNCC2.M, FG3DCIRCLE.M, ROT3Z.M, GR.M,  
#%                    FGRROT3.M, FRROT3.M, DRROT3.M.  
#% 
#% [x0n, an, rn, d, e, f, sigmah, conv, Vx0n, Van, urn, GNlog, a, ...  
#%   R0, R] = ls3dcircle(X, x0, a0, r0, tolp, tolg <, w >) 
#% --------------------------------------------------------------------- 

import numpy as np
from nlss11 import *
from gncc2 import *
from fg3dcircle import * 
from rot3z import * 
from gr import *  

def ls3dcircle(X=None, x0=None, a0=None, r0=None, tolp=None, tolg=None, w=None):
    """LS3DCIRCLE.py   Least-squares circle in three dimensions using Gauss-Newton. 
    Input     
    X    Array [x y z] where x = vector of x-coordinates,
         y = vector of y-coordinates and z = vector of z-coordinates.  
         Dimension: m x 3.
    x0   Estimate of the circle centre. 
         Dimension: 3 x 1. 
    a0   Estimate of the normal to the plane containing the circle.  
         Dimension: 3 x 1.  
    r0   Estimate of the circle radius.  
         Dimension: 1 x 1.  
    tolp Tolerance for test on step length.  
         Dimension: 1 x 1.  
    tolg Tolerance for test on gradient. 
         Dimension: 1 x 1.  
    Output   
    x0n  Estimate of the circle centre.  
         Dimension: 3 x 1.  
    an   Estimate of the normal direction  
         Dimension: 3 x 1. 
    rn   Estimate of the circle radius  
         Dimension: 1 x 1.
    d    Vector of distances from the points to the circle  
         Dimension: m x 1. 
    e    Vector of distances from the points to the plane  
         containing the circle 
         Dimension: m x 1.  
    f    Vector of distances from the points to the cylinder  
         containing the circle 
         Dimension: m x 1.  
    sigmah   Estimate of the standard deviation of the weighted
             residual errors.
             Dimension: 1 x 1.  
    conv     If conv = 1 the algorithm has converged,
             if conv = 0 the algorithm has not converged 
             and x0n, rn, d, e, f and sigmah are current estimates.  
             Dimension: 1 x 1.  
    Vx0n     Covariance matrix of circle centre.
             Dimension: 3 x 3.  
    Van      Covariance matrix of normal direction.  
             Dimension: 3 x 3.  
    urn      Uncertainty in circle radius.  
             Dimension: 1 x 1.  
    GNlog    Log of the Gauss-Newton iterations. 
             Rows 1 to niter contain 
             [iter, norm(f_iter), |step_iter|, |gradient_iter|]. 
             Row (niter + 1) contains  
             [conv, norm(d), 0, 0].
             Dimension: (niter + 1) x 4.  
    a        Optimisation parameters at the solution.
             Dimension: 6 x 1.  
    R0       Fixed rotation matrix. 
             Dimension: 3 x 3.  
    R        Upper-triangular factor of the Jacobian matrix 
             at the solution. 
             Dimension: 6 x 6. 
             """
    #check number of data points 
    if X==None: 
        raise ArithmeticError,'Require [x y z] array'
    else:
        if isinstance(X,np.matrix):
            pass
        else:
            x=np.mat(X);
        m = np.size(X, 0); 
        M1=np.mat(np.ones((m, 1)))
        if m < 6:
            raise ArithmeticError,"At least 6 data points required:"
    #if no weights are specified, use unit weights  
    if w==None:
        w = M1;
    #find the centroid of the data  
    xb = np.mat(np.mean(X,0).T);  
    #transform the data to close to standard position via a rotation  
    # followed by a translation 
    #R0 * a0 = [0 0 1]' 
    R0 = np.mat(rot3z(a0));
    #ת
    xb1 = R0 * xb;  
    #x
    x1 = R0 * x0;
    #n
    X1 = (np.mat(X) * R0.T);  
    #find xp, the point on axis nearest the centroid of the rotated data
    vz= np.mat([0,0,1]).T
    xp = x1 + (xb1[2,0] - x1[2,0]) * vz;  
    #translate data, mapping xp to the origin  
    X2 = X1 - M1 * xp.T;  
    x2 = x1 - xp;  
     
    ai = np.mat([0,0,0,0,0,r0]).T; 
    #tolp Tolerance for test on step length  
    #tolg Tolerance for test on gradient
    tol = np.mat([tolp,tolg]);  
      
    #Gauss-Newton algorithm to find estimate of roto-translation  
    #parameters that transform the data so that the best-fit circle  
    #is one in standard position 
    a, ef, R, GNlog= nlss11(ai, tol, 'fg3dcircle', X2, w);   
    e = ef[0:m];
    f = ef[m:2*m];
    d = np.sqrt(np.multiply(e,e) + np.multiply(f,f));  
    
    #inverse transformation to find circle centre and normal  
    #corresponding to original data
    rn = a[5]; 
    
    R3, DR1, DR2, DR3 = fgrrot3(np.concatenate((a[3:5],np.mat(0))));
    #an :axis  
    an = R0.T* R3.T * np.mat("0 ;0 ;1");
    #x0n:centre     
    x0n = R0.T* (xp + np.concatenate((a[0],a[1],a[2]),0)); 
    #GNlog    Log of the Gauss-Newton iterations.
    nGN = np.size(GNlog, 1);  
    conv = GNlog[nGN, 1];  
    if conv == 0 :
        raise ArithmeticError,'*** Gauss-Newton algorithm has not converged ***' 
     
    #calculate statistics  
    dof = 2 * m - 6;  
    sigmah = np.linalg.norm(d)/np.sqrt(dof);  
    G = np.mat(np.zeros((7, 6)));  
    G[0:3, 0] = R0.T * np.mat("1 0 0").T;  
    G[0:3, 1] = R0.T * np.mat("0 1 0").T;  
    G[0:3, 2] = R0.T * np.mat("0 0 1").T;  
    G[3:6, 3] = R0.T * DR1.T * np.mat("0 0 1").T;  
    G[3:6, 4] = R0.T * DR2.T * np.mat("0 0 1").T;  
    G[6, 5]   = 1; 
    #slove R' * Gt = sigmah * G'  
    Gt =np.linalg.solve(R.T,(sigmah * G.T));  
    Va = Gt.T * Gt; 
    #covariance matrix for x0n   
    Vx0n = Va[0:3, 0:3]; 
    #covariance matrix for an   
    Van = Va[3:6, 3:6]; 
    #uncertainty in rn  
    urn = np.sqrt(Va[6, 6]);
    
    #output x0n ,an,rn,d,e,f,sigmah,conv,Vx0n,Van,urn,GNlog,a,R0,R       
    return x0n, an, rn, d, e, f, sigmah, conv, Vx0n, Van, urn,GNlog, a, R0, R
if __name__=="__main__":
    print ls3dcircle.__doc__
    from lsplane import *
    from __log import *
    name=['x0', 'a', 'd', 'normd','x0n', 'an', 'rn', 'd', 'e', 'f', 'sigmah', 'conv', 'Vx0n', 'Van', 'urn','GNlog', 'a', 'R0', 'R']
    msg=''
    c11=np.mat([[182.86779073979     ,2325.87529134155     ,-294.40450739032 ],   
                [160.70160950720     ,2406.12086429678    , -294.58945521292] ,
                [213.76616244049    , 2504.70831446148    , -294.68201477650],
                [421.01919061425    , 2377.73490680325     ,-294.38556654694],
                [365.44126911669    , 2290.78641632270     ,-294.41985361918 ],
                [271.79572682184    , 2269.46402226983    , -294.34976761283 ],
                [222.96860277392    , 2287.43892234611     ,-294.44470758825 ]])
    
    res=lsplane(c11)
    for i in range(4):
        msg=msg+name[i]+":\n"+str(res[i])+"\n"
    res=ls3dcircle(c11, res[0], res[1], 54, 0.001, 0.001)
    
    msg=msg+"-------------------------------------------------\n"
    for i in range(4,len(res)+1):
        msg=msg+name[i]+":\n"+str(res[i-4])+"\n\n"
    msg=msg+"**************************************************\n"
    dlogger=initlog()
    dlogger.info(msg)
#ans
#x0n =
#
#  1.0e+003 *
#
#    0.2916
#    2.3992
#   -0.2945
#
#
#an =
#
#   -0.0004
#    0.0011
#    1.0000
#
#
#rn =
#
#  131.1465
#
#
#d =
#
#    0.0664
#    0.0347
#    0.0281
#    0.0476
#    0.0704
#    0.0568
#    0.0375
#
#
#e =
#
#    0.0016
#   -0.0298
#    0.0205
#    0.0008
#   -0.0305
#    0.0530
#   -0.0156
#
#
#f =
#
#    0.0664
#   -0.0177
#   -0.0192
#    0.0476
#   -0.0635
#    0.0204
#   -0.0341
#
#
#sigmah =
#
#    0.0478
#
#
#conv =
#
#     1
#
#
#Vx0n =
#
#  1.0e-003 *
#
#    0.7346    0.1964   -0.0001
#    0.1964    0.9803   -0.0004
#   -0.0001   -0.0004    0.5215
#
#
#Van =
#
#  1.0e-007 *
#
#    0.4272    0.1141    0.0000
#    0.1141    0.5698   -0.0006
#    0.0000   -0.0006    0.0000
#
#
#urn =
#
#    0.0228
#
#
#bGNlog =
#
#    1.0000  181.5656  184.2062  317.7424
#    2.0000   16.8206   19.8101   32.0230
#    3.0000    0.1546    0.0941    0.1300
#    4.0000    0.1351    0.0000    0.0000
#    1.0000    0.1351         0         0
#
#
#a =
#
#   28.9816
#   47.4224
#    0.0000
#   -0.0000
#    0.0000
#  131.1465
#
#
#R0 =
#
#    1.0000    0.0000    0.0004
#         0    1.0000   -0.0011
#   -0.0004    0.0011    1.0000
#
#
#R =
#
#   -1.9038    0.0516    0.0000   -0.0266    0.0203    0.8125
#         0    1.8365   -0.0000   -0.0203   -0.0281   -1.4011
#         0         0    2.6458  125.4678  -76.6782   -0.0000
#         0         0         0 -205.7199  -54.9671    0.0000
#         0         0         0         0  231.1537   -0.0002
#         0         0         0         0         0    2.0921

  
