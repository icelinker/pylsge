#tested beta
#  function [dR1, dR2, dR3] = drrot3(R1, R2, R3) 
#% -------------------------------------------------------------------------- 
#% DRROT3.M   Calculate the derivatives of plane rotations -  
#%            use right-handed rotation matrices. 
#% 
#% Version 1.0     
#% Last amended   I M Smith 27 May 2002.  
#% Created        I M Smith 08 Mar 2002 
#% -------------------------------------------------------------------------- 
#% Input  
#% R1       Plane rotation of the form 
#%          [1 0 0; 0 c1 -s1; 0 s1 c1].  
#%          Dimension: 3 x 3.  
#%  
#%  
#% 
#% Output  
#% dR1      Derivative of R1 with respect to rotation angle.  
#%          [0 0 0; 0 -s1 -c1; 0 c1 -s1].  
#%          Dimension: 3 x 3.  
#%  
#%  
#%  
#% [dR1 <, dR2, dR3 >] = drrot3(R1 <, R2, R3 >) 
#% -------------------------------------------------------------------------- 
# 
#  if nargin > 0 
#    dR1 = [0 0 0; 0 -R1(3, 2) -R1(2, 2); 0 R1(2, 2) -R1(3, 2)]; 
#  end % if nargin  
#% 
#  if nargin > 1 
#    dR2 = [-R2(1, 3) 0 R2(1, 1); 0 0 0; -R2(1, 1) 0 -R2(1, 3)]; 
#  end % if nargin  
#% 
#  if nargin > 2 
#    dR3 = [-R3(2, 1) -R3(1, 1) 0; R3(1, 1) -R3(2, 1) 0; 0 0 0]; 
#  end % if nargin  
#% -------------------------------------------------------------------------- 
#% End of DRROT3.M.
import numpy as np

def drrot3(R1=None, R2=None, R3=None):
    """ DRROT3.M   Calculate the derivatives of plane rotations use right-handed rotation matrices. 
    Input:  
    R1    Plane rotation of the form  
    [1 0 0; 0 c1 -s1; 0 s1 c1].  
    Dimension: 3 x 3.  
    Output:  
    dR1    Derivative of R1 with respect to rotation angle.  
    [0 0 0; 0 -s1 -c1; 0 c1 -s1].  
    Dimension: 3 x 3.  
    [dR1 <, dR2, dR3 >] = drrot3(R1 <, R2, R3 >) """
    if isinstance(R1,np.matrix):
        R1=R1;
    else:
        R1=np.mat(R1);
    if R1!=None:    
        dR1 =np.mat([[0,0,0],
                     [ 0,-R1[2, 1],-R1[1,1]],
                     [ 0, R1[1, 1],-R1[2, 1]]]); 
    if R2!=None:
        dR2 = np.mat([[-R2[0, 2],0,R2[0, 0]],
               [ 0,0,0],
               [ -R2[0, 0], 0 ,-R2[0, 2]]]);    
    if R3!=None: 
        dR3 =np.mat([[-R3[1, 0],-R3[0, 0],0],
                     [ R3[0, 0],-R3[1, 0],0],
                     [ 0,0,0]]);
    if R1!=None and R2!=None and R3!=None :
        return  dR1,dR2,dR3
    elif R1!=None and R2!=None:
        return  dR1,dR2,None
    elif R1!=None:
        return  dR1,None,None
if __name__=="__main__":
    print drrot3.__doc__ 
    x=np.mat("1 3 2;4 5 6;7 8 9")
    dR1,dR2, dR3 = drrot3(x,x,x)
    print "dR1=",dR1,"\ndR2=",dR2, "\ndR3=",dR3
    dR1,dR2, dR3 = drrot3(x,x)
    print "dR1=",dR1,"\ndR2=",dR2, "\ndR3=",dR3
    dR1,dR2, dR3 = drrot3(x)
    print "dR1=",dR1,"\ndR2=",dR2, "\ndR3=",dR3






