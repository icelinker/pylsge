#  function [c, iconv, sp, sg] = gncc2(f0, f1, p, g, scale, tolr, scalef) 
#% -------------------------------------------------------------------------- 
#% GNCC2   Gauss-Newton convergence conditions. 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% -------------------------------------------------------------------------- 
#% Input  
#% f0       Norm of residuals at old point: norm(f(a))  
#%          Dimension: 1 x 1. 
#%  
#% f1       Norm of residuals at new point: norm(f(a + p)) 
#%          Dimension: 1 x 1. 
#% 
#% p        Step 
#%          Dimension: n x 1. 
#%  
#% g        Gradient 
#%          Dimension: n x 1. 
#%  
#% scale    Scale for columns of Jacobian matrix. 
#%          Dimension: n x 1. 
#%  
#% tolr     Relative tolerance. 
#%          Dimension: 1 x 1. 
#%  
#% scalef   Scale for function values.  
#%          Dimension: 1 x 1. 
#% 
#% Output  
#% c        Convergence indices in the form of ratio of value over 
#%          scaled tolerance. 
#%          c(1) size of step 
#%          c(2) change in function value 
#%          c(3) size of gradient 
#%          c(4) sum of squares near zero 
#%          c(5) gradient near zero 
#%          Dimension: 1 x 5. 
#%  
#% iconv    = 1 if convergence criteria have been met, i.e.,  
#%          C(1), C(2), C(3) < 1 or  
#%                      C(4) < 1 or  
#%                      C(5) < 1. 
#%          = 0 otherwise. 
#%          Dimension: 1 x 1. 
#% 
#% sp       Scaled size of the step. 
#%          Dimension: 1 x 1. 
#%  
#% sg       Scaled size of the gradient. 
#%          Dimension: 1 x 1. 
#% 
#% [c, iconv, sp, sg] = gncc2(f0, f1, p, g, scale, tolr, scalef) 
#% --------------------------------------------------------------------------  

import numpy as np
def gncc2(f0=None, f1=None, p=None, g=None, scale=None, tolr=None, scalef=None):
    """
     GNCC2.py   Gauss-Newton convergence conditions. 
     Input  
        f0    Norm of residuals at old point: norm(f(a))  
              Dimension: 1 x 1. 
        f1    Norm of residuals at new point: norm(f(a + p)) 
              Dimension: 1 x 1. 
        p     Step 
              Dimension: n x 1.
        g     Gradient 
              Dimension: n x 1. 
        scale    Scale for columns of Jacobian matrix. 
                  Dimension: n x 1. 
        tolr     Relative tolerance. 
        Dimension: 1 x 1. 
        scalef   Scale for function values.  
          Dimension: 1 x 1. 
        Output  
         c        Convergence indices in the form of ratio of value over 
                  scaled tolerance. 
                  c(1) size of step 
                  c(2) change in function value 
                  c(3) size of gradient 
                  c(4) sum of squares near zero 
                  c(5) gradient near zero 
                  Dimension: 1 x 5. 
  
         iconv    = 1 if convergence criteria have been met, i.e.,  
                     C(1), C(2), C(3) < 1 or  
                    C(4) < 1 or  
                      C(5) < 1. 
          = 0 otherwise. 
          Dimension: 1 x 1. 
 
         sp       Scaled size of the step. 
              Dimension: 1 x 1. 
  
         sg       Scaled size of the gradient. 
              Dimension: 1 x 1. 
 
         [c, iconv, sp, sg] = gncc2(f0, f1, p, g, scale, tolr, scalef) """
    iconv = 0; 
    eps=np.finfo(np.double).eps
    c=np.zeros((5,1))
    #scale Scale for columns of Jacobian matrix.
    
    sp = np.absolute(np.multiply(p,scale)).max(); 
    sg = np.absolute(g/scale).max(); 
    c[0] = sp/(scalef * tolr**(0.7)); 
    #print tolr ,scalef
    delf = f0 - f1; 
    #print delf
    c[1]= np.absolute(delf)/(tolr * scalef); 
    
    d3 = (tolr**(0.7)) * (scalef); 
    d4 = scalef * (eps**(0.7)); 
    #print d4
    d5 = (eps**(0.7)) * (scalef); 
    #print d5
    
    c[2] = sg/d3; 
    c[3] = f1/d4; 
    c[4] = sg/d5; 
    
    if c[0] < 1.0 and  c[1] < 1.0 and c[2] < 1.0 :
        iconv = 1; 
    elif (c[3]) < 1.0 : 
        iconv = 1; 
    elif (c[4] < 1.0): 
        iconv = 1;
        
    return c, iconv, sp, sg 
if __name__=="__main__":
    print gncc2.__doc__
    a=np.mat("1.1 2.2 1.3 4.4 1.2").T
    F0=181.56558632164223
    F1=16.820616190420548
    p=np.mat("27.24694215746634;\
            44.344806248245305;\
            -0.000000000000016;\
            0.000001157097032;\
            -0.000007741612598;\
            69.62339308557026")
    g=np.mat("-76.2245613814803;\
            -116.41405816408438;\
            0.000000000000226;\
            0.181841871004714;\
            3.150320092073721;\
            -840.6674088836138")
    scale=np.mat("1.943645896458876;\
            1.795060062833155;\
            2.645751311064591;\
            205.71989478451871;\
            237.59929946421013;\
            2.645751311064591")
    tol=np.mat("0.001    0.001")
    #c:23190184.546423994    164744970.1312217    40001399.96093918    1525224396859369.2    28811577205613188
    #conv:0
    #sp:184.20618353691287
    #sg:317.7424141747276
    c,conv,sp,sg=gncc2(F0, F1, p, g, scale, tol[0,0], tol[0,1])
    print "c=",c,"\nconv=",conv,"\nsp=",sp,"\nsg=",sg





