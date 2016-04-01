#tested beta
# function [R, R1, R2, R3] = frrot3(theta, U0) 
#% -------------------------------------------------------------------------- 
#% FRROT3.M   Form rotation matrix R = R3*R2*R1*U0. - use right-handed 
#%            rotation matrices. 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% -------------------------------------------------------------------------- 
#% Input  
#% theta    Array of plane rotation angles (t1, t2, t3). 
#%          Dimension: 3 x 1.  
#%  
#%  
#% 
#% Output  
#% R        Rotation matrix.  
#%          Dimension: 3 x 3.  
#%  
#% R1       Plane rotation [1 0 0; 0 c1 -s1; 0 s1 c1]. 
#%	       Dimension: 3 x 3.  
#%  
#%  
#% 
#% [R, R1 <, R2, R3 >] = frrot3(theta <, U0 >) 
#% -------------------------------------------------------------------------- 
# 
#  ct = cos(theta);  
#  st = sin(theta); 
#% 
#  if length(theta) > 0 
#    R1 = [ 1 0 0; 0 ct(1) -st(1); 0 st(1) ct(1)]; 
#    R = R1; 
#  end %if 
#% 
#  if length(theta) > 1 
#    R2 = [ ct(2) 0 st(2); 0 1 0; -st(2) 0 ct(2)]; 
#    R = R2*R; 
#  end % if 
#% 
#  if length(theta) > 2 
#    R3 = [ ct(3) -st(3) 0; st(3) ct(3) 0; 0 0 1]; 
#    R = R3*R; 
#   end % if 
#% 
#  if nargin > 1 
#    R = R*U0; 
#  end % if 
#% -------------------------------------------------------------------------- 
#% End of FRROT3.M.
import numpy as np

def frrot3(theta=None, U0=None): 
    """FRROT3.py   Form rotation matrix R = R3*R2*R1*U0. - use right-handed rotation matrices. 
    Input  
    theta Array of plane rotation angles (t1, t2, t3). 
          Dimension: 3 x 1. 
    U0    Array of plane rotation. 
          Dimension: 3 x 1.  
    Output  
    R     Rotation matrix.  
          Dimension: 3 x 3.  
    R1    Plane rotation [1 0 0; 0 c1 -s1; 0 s1 c1]. 
          Dimension: 3 x 3.  
    [R, R1 <, R2, R3 >] = frrot3(theta <, U0 >) """
    R, R1,R2, R3=None,None,None,None
    if theta==None:
        raise  TypeError,"theta type error"
    elif isinstance(theta,np.ndarray)==False:
        print "covert theta type "
        try :
            theta=np.asmatrix(theta)
        except:
            print "theta type covert error"
            return None
    ct = np.cos(theta); 
    #print ct 
    st = np.sin(theta);
    #print st 
    if theta.size > 0:
        R1 = np.mat([[1,0,0], 
                     [0,ct[0,0],-st[0,0]],
                    [0,st[0,0],ct[0,0]]]); 
        R = R1; 
        #print R
    
    if theta.size > 1 :
        R2 = np.mat([[ct[1,0], 0, st[1,0]],
                    [0 ,1, 0],
                    [-st[1,0] ,0 ,ct[1,0]]]);
        R = R2*R;
        #print R 
    if theta.size > 2 :
        R3 = np.mat([[ ct[2,0] ,-st[2,0], 0],
                     [st[2,0], ct[2,0], 0],
                     [0 ,0 ,1]]);
        R = R3*R; 
    if U0!=None:
        if isinstance(U0,np.ndarray)==False:
            try :
                U0=np.asmatrix(U0)
            except:
                print "U0 type covert error"
                return None
        R = R*U0;
    else:
        R = R;
    return R, R1,R2, R3
    
if __name__=="__main__":
    import math
    print frrot3.__doc__
    theta=np.mat([[ 1.40048515],
                  [ 2.10329182],
                  [ 1.12845222]]);
    print theta
    U0=np.mat(np.eye(3))
    R, R1,R2, R3 = frrot3(theta ,U0 )
    print "R=",R,"\nR1=",R1,"\nR2=",R2,"\nR3=",R3






