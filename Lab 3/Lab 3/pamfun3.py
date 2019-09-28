# File: pamfun.py 
# Functions for pulse amplitude modulation (PAM) import numpy as np import comsig
import numpy as np
import comsig
def pam10(sps, ptype, pparms=[]): 
     
    # ***** Set up PAM pulse p(t) ***** 
    ptype = ptype.lower() # Convert ptype to lowercase
                          # Set left/right limits for p(t) 
    if (ptype=="rect"): 
        kL = -0.5; kR = -kL 
    else: 
        kL = -1.0; kR = -kL # Default left/right limits 
    ixpL = int(np.ceil(sps*kL) # Left index for p(t) time axis 
    ixpR = int(np.ceil(sps*kR) # Right index for p(t) time axis 
    ttp = np.arange(ixpL,ixpR) # Time axis for p(t) 
    pt = np.zeros(len(ttp)) # Initialize pulse p(t) 
    if (ptype=="rect"): # Rectangular p(t) 
        pt = np.ones(ttp)
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
    else:
        print("ptype ’%s’ is not recognized" % ptype)
    # ***** Filter with h(t) = p(t) ***** 
    return(pt)