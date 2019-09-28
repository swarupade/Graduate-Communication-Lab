# this module will be imported in the into your flowgraph
# File: pamfun.py 
# Functions for pulse amplitude modulation (PAM) import numpy as np import comsig
import numpy as np

def pam10(sps, ptype, pparms=[]): 
     
    # ***** Set up PAM pulse p(t) ***** 
    ptype = ptype.lower() # Convert ptype to lowercase
                          # Set left/right limits for p(t) 
    if (ptype=="rect"): 
        kL = -0.5; kR = -kL 
    else: 
        kL = -1.0; kR = -kL # Default left/right limits 
    
    ixpL = int(np.ceil(sps*kL)) # Left index for p(t) time axis 
    ixpR = int(np.ceil(sps*kR)) # Right index for p(t) time axis 
    ttp = np.arange(ixpL,ixpR) # Time axis for p(t) 
    pt = np.zeros(len(ttp)) # Initialize pulse p(t) 
    if (ptype=="rect"): # Rectangular p(t) 
        pt = np.ones(len(ttp))
    elif (ptype == "tri"):
	    pt = np.array([(1+i/sps) if i/sps+1 < 1.0 else (1-i/sps) for i in list(ttp)]) 
    elif (ptype=="sinc"):
        k=pparms[0]
        kL = -1.0*k
        kR = -kL
        ixpL = int(np.ceil(sps*kL)) 		# Left index for p(t) time axis
        ixpR = int(np.ceil(sps*kR))		# Right index for p(t) time axis
        ttp = np.arange(ixpL,ixpR)	 # Time axis for p(t)
        pt = np.zeros(len(ttp))
        beta=pparms[1]
        pt = np.array([np.sin(np.pi*t/sps)/(np.pi*t/sps) if t!= 0  else  1.0 for t in list(ttp)])
        pt=pt*np.kaiser(len(pt), beta)
    elif (ptype=="man"):
        ix = np.where(np.logical_and(ttp>=kL*sps, ttp<=kR*sps))[0]
        pt[ix] = np.ones(len(ix))
        pt[ix] = np.concatenate([-1*np.ones(int(len(ix)/2)),np.ones(len(pt)-int(len(pt)/2))])
        #for i in ix:
            #print(pt[i])
    
    elif (ptype=='rcf'):
        k= pparms[0]
        alpha=pparms[1]
        kL = -k; kR = -kL
        ixpL = np.ceil(sps*kL)		# Left index for p(t) time axis
        ixpR = np.ceil(sps*kR)		# Right index for p(t) time axis
        ttp = np.arange(ixpL,ixpR)/float(sps)	 # Time axis for p(t)
        pt = []
        
        for t in ttp:
            num = (np.sin(np.pi*t))*(np.cos(np.pi*alpha*t))
            denum =(np.pi*t)*(1-np.power(2*alpha*t,2))
            rcf_t = np.divide(num,(float(denum)))
            if denum == 0.0:
                rcf_t = 1
            pt = np.concatenate([pt,[rcf_t]])
    else:
        print("ptype is not recognized")
    # ***** Filter with h(t) = p(t) ***** 
    return(pt)
