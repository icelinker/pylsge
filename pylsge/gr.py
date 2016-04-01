#tested beta
# function [U, c, s] = gr(x, y) 
#% -------------------------------------------------------------------------- 
#% GR.M   Form Givens plane rotation.  
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% -------------------------------------------------------------------------- 
#% Input  
#% x        Scalar. 
#%          Dimension: 1 x 1.  
#%  
#% y        Scalar. 
#%          Dimension: 1 x 1.  
#%  
#% Output  
#% U        Rotation matrix [c s; -s c], with U * [x y]' = [z 0]'.  
#%          Dimension: 2 x 2.  
#%  
#% c        Cosine of the rotation angle.  
#%          Dimension: 1 x 1.  
#%  
#% s        Sine of the rotation angle.  
#%          Dimension: 1 x 1.  
#%  
#% [U, c, s] = gr(x, y) 
#% -------------------------------------------------------------------------- 
# 
#% form sine and cosine: s and c 
#  if y == 0 
#    c = 1; 
#    s = 0; 
#  elseif  abs(y) >= abs(x) 
#    t = x/y; 
#    s = 1/sqrt(1 + t*t); 
#    c = t*s; 
#  else 
#    t = y/x; 
#    c = 1/sqrt(1 + t*t); 
#    s = t*c; 
#  end % if 
#% 
#% assign U 
#  U = [c  s; -s  c]; 
#% -------------------------------------------------------------------------- 
#% End of GR.M.
import numpy as np
def gr(x=None,y=None):
    """    gr.py   Form Givens plane rotation.
    [U, c, s] = gr(x, y) 
    Input  
    x        Scalar. 
          Dimension: 1 x 1.  
    y        Scalar. 
          Dimension: 1 x 1.  
    Output  
    U        Rotation matrix [c s; -s c], with U * np.array[x y].transpose() = [r 0]'.  
          Dimension: 2 x 2.  
    
    c        Cosine of the rotation angle.  
          Dimension: 1 x 1.  
          
    s        Sine of the rotation angle.  
          Dimension: 1 x 1. """ 
    #form sine and cosine: s and c 
    if x!=None and y!=None:
        if y == 0.0 and x==0.0: 
            c = 1.0; 
            s = 0.0; 
        r=np.sqrt(x*x+y*y)
        s =y/r
        c =x/r;
        U=np.mat([[c,s],[-s,c]])
        return  U,c,s;
    else:
        raise  TypeError,"xy type error"
if __name__=="__main__":
    print gr.__doc__
    r,c,s=gr(1.0,1.0)
    print r
    r,c,s=gr(1.0,1.0)
    g= r*np.mat([1.0,1.0]).T
    print g
    r,s,c=gr(1.0,-10.0)
    print r[0]
    