#coding:utf-8
# function [a, f, R, GNlog] = nlss11(ai,tol,fguser,p1,p2,p3,p4,p5,p6, ...  
#            p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20) 
#% -------------------------------------------------------------------------- 
#% NLSS11.M   Nonlinear least squares solver: 
#%            Minimise f'*f where  
#%            [f, J] = fguser(a, p1, p2,...). 
#% 
#% Author   A B Forbes, CMSC, NPL. 
#% Amendment history 
#% v1.1  2002-05-03 IMS Prepared for MetroS.  
#% v1.1  2002-01-04 ABF Output a, f, GNlog set to latest value. 
#% v1.1  2001-07-09 ABF Optional input arguments increased. 
#% v1.1  2000-12-19 ABF Better convergence criteria. Column scaling. 
#% v1.0a 1999-07-16 ABF Statistics removed. 
#% v1.0a 1999-07-16 ABF Created. 
#% -------------------------------------------------------------------------- 
#% Input  
#% ai       Optimisation parameters, intial estimates. 
#%          Dimension: n x 1 
#% 
#% tol      Convergence tolerances [tolr tols]', where  
#%          tolr = relative tolerance, and  
#%          tols = scale for function values.  
#%          Dimension: 2 x 1.  
#% 
#% fguser   Module to calculate function and gradients with 
#%          signature 
#%          [f, J] = fguser(a, p1, p2,...). 
#% 
#% p1,...   Additional parameters to be passed to fguser  
#%          without change. 
#%          Note: The number of additional parameters supported is 20.  
#%            
#% Output  
#% a        Solution estimates of the optimisation parameters. 
#%          Dimension: n x 1. 
#% 
#% f        Functions evaluated at a. 
#%          Dimension: m x 1. 
#%          Constraint: m >= n. 
#%  
#% R        Triangular factor of the Jacobian matrix evaluated at a. 
#%          Dimension: n x n. 
#% 
#% GNlog    Log of the Gauss-Newton iterations.  
#%          Rows 1 to niter contain  
#%          [iter, norm(f_iter), |step_iter|, |gradient_iter|].  
#%          Row (niter + 1) contains  
#%          [conv, norm(d), 0, 0].  
#%          Dimension: (niter + 1) x 4.  
#%  
#% Modular structure: GNCC2.M.  
#%  
#% [a, f, R, GNlog] = nlss11(a, tol, 'fguser', p1, p2, p3, ...) 
#% -------------------------------------------------------------------------- 

from fg3dcircle import*
from gncc2 import *
def nlss11(ai=None, tol=None, spline_type='fguser', *args): 
    """NLSS11.py   Nonlinear least squares solver: 
            Minimise f'*f where [f, J] = fguser(a, p1, p2,...). 
        Input  
        ai       Optimisation parameters, intial estimates. 
                  Dimension: n x 1 
        tol      Convergence tolerances [tolr tols]', where  
                  tolr = relative tolerance, and  
                  tols = scale for function values.  
                  Dimension: 2 x 1.  
        fguser   Module to calculate function and gradients with signature 
                  [f, J] = fguser(a, p1, p2,...). 
        p1,...   Additional parameters to be passed to fguser  
                  without change. 
                  Note: The number of additional parameters supported is 20.  
        Output  
        a        Solution estimates of the optimisation parameters. 
                  Dimension: n x 1. 
        f        Functions evaluated at a. 
                  Dimension: m x 1. 
                  Constraint: m >= n. 
        R        Triangular factor of the Jacobian matrix evaluated at a. 
                  Dimension: n x n. 
        GNlog    Log of the Gauss-Newton iterations.  
                  Rows 1 to niter contain  
                  [iter, norm(f_iter), |step_iter|, |gradient_iter|].  
                  Row (niter + 1) contains  
                  [conv, norm(d), 0, 0].  
                  Dimension: (niter + 1) x 4.  
        Modular structure: GNCC2.py
         a, f, R, GNlog = nlss11(a, tol, 'fguser', p1, p2, p3, ...) """ 
    spline_gr={'fguser':None,
             'fg3dcircle':fg3dcircle,
             }    
    a0 = ai; 
    n = np.size(a0,0); 
    if n == 0: 
        raise ArithmeticError,"Empty vector of parameter estimates:"
    if tol!=None:
        tol=np.asarray(tol)
    # Set up call to fguser: 
#    f0, J = fguser(a0,p1,p2,....);
#    if spline_gr[spline_type] !=None:
#        pass
#    else:
#        raise ArithmeticError,"Unkown Spline_type"
#
#    callfg0 = [fguser,'(a0',]; 
#    for k = 1:nargin-3 :
#        callfg0 = [callfg0,',p',int2str(k),]; 
##  end % for k 
#    callfg0 = [callfg0,')'];      
##% 
#    f1,_= fguser(a1,p1,p2,....); 
#    callfg1 = [fguser,'(a1',]; 
#    for k = 1:nargin-3 :
#        callfg1 = [callfg1,',p',int2str(k),]; 
##  end % for k   
#    callfg1 = [callfg1,')'];      
#% 
    fguser=spline_gr['fg3dcircle']
    mxiter = int(100+np.ceil(np.sqrt(n))); 
    conv = 0; 
    niter = 0; 
    eta = 0.01; 
    GNlog = []; 
    #% G-N iterations 
    while niter < mxiter and conv == 0 :
        f0, J = fguser(a0,*args);
        if niter == 0: #scale by norm of columns of J :
            mJ,nJ = J.shape; 
            scale = np.zeros((nJ,1)); 
            for j in range(nJ): 
                scale[j] = np.linalg.norm(J[:,j]); 
            #end 
        #end % if nite
        m = np.size(f0,0); 
        #% Check on m, n. 
        if niter == 0 and m < n :
            raise ArithmeticError ,"Number of observation less than number of parameters" 
        #endif
          
        #% Calculate update step and gradient. 
        F0 = np.linalg.norm(f0); 
        Rj= np.concatenate((J, f0),1);
        Ra = np.triu(np.linalg.qr(Rj,mode = 'economic')); 
        R = np.mat(Ra[0:nJ,0:nJ]); 
        q = np.mat(Ra[0:nJ,nJ]).T; 
        p = np.linalg.solve(-R,q); 
        g = 2*R.T*q; 
        G0 = np.asscalar(g.T*p); 
        a1 = a0 + p; 
        niter = niter + 1; 
        #% 
        #% Check on convergence. 
        #[f1] = eval(callfg1); 
        f1, J = fguser(a0,*args);
        F1 = np.linalg.norm(f1); 
        c, conv, sp, sg = gncc2(F0, F1, p, g, scale, tol[0,0], tol[0,1]); 
        #testlog.error( c, conv, sp, sg)
        if conv != 1 :
            #% ...otherwise check on reduction of sum of squares.
            #% Evaluate f at a1. 
            rho = (F1 - F0)*(F1 + F0)/(G0); 
            if rho < eta :
                tmin = np.maximum(0.001,\
                                  1/(2*(1-rho))); 
                a0 = a0 + tmin*p; 
            else: 
                a0 = a0 + p; 
            #end % if rho 
        #end % if conv  
        GNlog.append([niter, F0, sp, sg]);
    #  end % while niter  
    a = a0+p; 
    f = f1; 
    GNlog.append([conv, F1, 0, 0]); 
    return a, f, R, np.mat(GNlog)
    
if __name__=="__main__":
    print nlss11.__doc__
    J=np.mat("0.951286979964784    0.308306798091578    0    0.020479239790891    -0.06318911646854    -1;\
    0.882298966025384    -0.470689424727749    0    0.008344755411305    0.01564209579467    -1;\
    0.304397702133556    -0.952545032497578    0    0.018272619293798    0.005839244482132    -1;\
    -0.986787652935496    -0.162018912519668    0    -0.007712309728096    0.046972368206556    -1;\
    -0.860167393151313    0.510011819234882    0    -0.03236419346161    -0.054584272111734    -1;\
    -0.110470744211805    0.993879376319526    0    0.020289233466095    0.00225516976596    -1;\
    0.525225117929332    0.850963322062778    0    -0.028993065800175    0.017894879848778    -1;\
    0    0    -1    -25.857441170012407    79.78366767291251    0;\
    0    0    -1    54.388290807379235    101.94988504573811    0;\
    0    0    -1    152.9757823581317    48.88532829760638    0;\
    0    0    -1    26.002118728390997    -158.36774431024708    0;\
    0    0    -1    -60.94627570031207    -102.78977293127235    0;\
    0    0    -1    -82.26873595533516    -9.144257042545291    0;\
    0    0    -1    -64.29373906824367    39.68289326780737    0")
    f0=np.mat("29.86918916504675 ;61.5502714743123;106.5969031795048;106.48816970817859;\
    65.49973197041481;28.775372862638278;21.554066081711284;0.06642487261928;-0.017728793070148;\
    -0.019182945341583;0.04760129301053;-0.063457732234838;0.020414180985654;-0.034070875969007")
    Rj= np.concatenate((J, f0),1);
    Ra = np.triu(np.linalg.qr(Rj,mode = 'economic'));
    ai=np.mat("0;0;0;0;0;54")
    tol=np.mat("0.001  0.001")
    X2=np.mat("-79.78366767291251    -25.857441170012407    0.06642487261928;\
            -101.94988504573811    54.388290807379235    -0.017728793070148;\
            -48.88532829760638    152.9757823581317    -0.019182945341583;\
            158.36774431024708    26.002118728390997    0.04760129301053;\
            102.78977293127235    -60.94627570031207    -0.063457732234838;\
            9.144257042545291    -82.26873595533516    0.020414180985654;\
            -39.68289326780737    -64.29373906824367    -0.034070875969007")
    a,ef,R,GNlog=nlss11(ai,tol,'fg3dcircle',X2)
    print "a=",a,"\nef=",ef,"\nR=",R,"\nGNlog=",GNlog
#a:[28.981619128152957
#47.42237593383242
#0.000004384459277
#-0.0000000852615
#0.000000011771474
#131.14648073390597]

#ef:[0.001600482293526
#-0.02980461947439
#0.020523916154644
#0.000752589389634
#-0.030470171160857
#0.053006449173608
#-0.015608646375824
#0.066428016434899
#-0.017732230196964
#-0.019195412833883
#0.047597211809092
#-0.063453745851077
#0.020421087700129
#-0.034064927062196]

##R:-1.903849952408756    0.051585367077737    0.000000025253554    -0.026576563822651    0.020267314998876    0.81251487571466
#0    1.836489670163342    -0.000000169611155    -0.020280895435765    -0.028113359940918    -1.401098842345307
#0    0    2.645751311064579    125.46780320653346    -76.67815511637673    -0.00000018602004
#0    0    0    -205.71989813525653    -54.96712547025054    0.000033046605773
#0    0    0    0    231.15371190179152    -0.000233848200741
#0    0    0    0    0    2.092066336171952


#GNLog:1    181.56558632164223    184.20618353691287    317.7424141747276
#2    16.820616190420548    19.810063564984905    32.02297200633283
#3    0.154608556557191    0.094147872631575    0.130034889093287
#4    0.135127010718464    0.000007571221536    0.000007095796329
#1    0.135127010546328    0    0



    


    








