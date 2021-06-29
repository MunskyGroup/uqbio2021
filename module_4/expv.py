import numpy as np
import scipy.sparse.linalg as spl

def expv(t,A,v,tol=1.0e-7,m=30):
    '''
    a python version of expv from 
    roger sidje's expokit. A should be sparse.
    '''
    n,n = A.shape 
    if n<m:
        m = n
    anorm = spl.norm(A,np.inf) 

    # set some default parameters
    mxrej = 10.0;  btol  = 1.0e-7; 
    gamma = 0.9; delta = 1.2; 
    mb    = m; t_out   = abs(t);
    nstep = 0.0; t_new   = 0.0;
    t_now = 0.0; s_error = 0.0;
    eps = np.finfo(float).eps 
    rndoff= anorm*eps; 

    k1 = 2; xm = 1.0/m; normv = np.linalg.norm(v); beta = normv;
    fact = (((m+1)/np.exp(1))**(m+1))*np.sqrt(2*np.pi*(m+1));
    t_new = (1.0/anorm)*((fact*tol)/(4.0*beta*anorm))**xm;
    s = 10**(np.floor(np.log10(t_new))-1); t_new = np.ceil(float(t_new)/s)*s; 
    sgn = np.sign(t); nstep = 0;   
    w = v
    hump  = normv
    while t_now < t_out:
        nstep += 1
        t_step = np.min([t_out-t_now,t_new])
        V = np.zeros((n,m+1))
        H = np.zeros((m+2,m+2))
        V[:,0] = np.ravel((1.0/beta)*w)
        for j in range(m):
            p = A.dot(V[:,j])
            for i in range(j+1): 
                H[i,j] = np.dot(V[:,i].T,p)
                p = p-H[i,j]*V[:,i]
            s = np.linalg.norm(p)
            if s<btol: 
                k1 = 0.0
                mb = j+1
                t_step = t_out - t_now
                break                
            H[j+1,j] = s
            V[:,j+1] = (1/s)*p
        if k1 != 0:
            H[m+1,m] = 1
            avnorm = np.linalg.norm(A.dot(V[:,m])) # check index.
        ireject = 0
        while ireject<=mxrej:
            mx = int(mb+k1)
            F = spl.expm(sgn*t_step*H[:mx,:mx])
            if k1 == 0:
                err_loc = btol
                break
            else:
                phi1 = np.abs(beta*F[m,1])
                phi2 = np.abs(beta*F[m+1,1] * avnorm)
                if phi1 > 10*phi2:
                    err_loc = phi2
                    xm = 1.0/float(m)
                elif phi1>phi2:
                    err_loc = (phi1*phi2)/(phi1-phi2)
                    xm = 1.0/float(m)
                else:
                    err_loc = phi1
                    xm = 1.0 / float(m-1)
                
            if err_loc <= delta * t_step * tol:
                break
            else:
                t_step = gamma * t_step * (t_step * tol/err_loc)**xm
                s = 10**(np.floor(np.log10(t_step))-1)
                t_step = np.ceil(t_step/s)*s
                if ireject == mxrej:
                    print('The requested tolerance is too high.')
                ireject += 1
      
        mx = int(mb) + int(np.max( [0,k1-1]))                    
        w = np.dot(V[:,:mx],beta*F[:mx,0])
        beta = np.linalg.norm(w)
        hump = np.max([hump,beta])
        t_now += t_step
        t_new = gamma * t_step * (t_step*tol/float(err_loc))**xm
        s = 10**(np.floor(np.log10(t_new))-1)
        t_new = np.ceil(t_new/float(s))*s
        err_loc = np.max([err_loc,rndoff])
        s_error = s_error + err_loc
    err = s_error
    hump = hump / float(normv)
    return w, err, hump
        
