# %%writefile pam11.py
# File: pam11.py
# Functions for pulse amplitude modulation (PAM)
from pylab import *
import comsig 
import matplotlib as plt
def pam11(sig_an, FB, Fs, ptype, pparms=[]):
    """
    Pulse amplitude modulation: a_n -> s(t), -TB/2<=t<(N-1/2)*TB,
    V1.1 for ’man’, ’rcf’, ’rect’, ’sinc’, and ’tri’ pulse types.
    >>>>> sig_st = pam11(sig_an, Fs, ptype, pparms) <<<<<
    where  sig_an: sequence from class sigSequ
        sig_an.signal():  N-symbol DT input sequence a_n, 0 <= n < N
        sig_an.get_FB():  Baud rate of a_n, TB=1/FB
        Fs:    sampling rate of s(t)
        ptype: pulse type from list
        (’man’,’rcf’,’rect’,’sinc’,’tri’)
        pparms not used for ’man’,’rect’,’tri’
        pparms = [k, alpha] for ’rcf’
        pparms = [k, beta]  for ’sinc’
        k:     "tail" truncation parameter for ’rcf’,’sinc’
        (truncates p(t) to -k*TB <= t < k*TB)
        alpha: Rolloff parameter for ’rcf’, 0<=alpha<=1
        beta:  Kaiser window parameter for ’sinc’
        sig_st: waveform from class sigWave
        sig_st.timeAxis():  time axis for s(t), starts at -TB/2
        sig_st.signal():    CT output signal s(t), -TB/2<=t<(N-1/2)*TB,
        with sampling rate Fs
    """

    # ***** Set variables and manage data formatting *****
    
    an = sig_an.signal()
    N = len(an)    		# Number of data symbols
    TB = 1/float(FB)    	# Time per symbol
    ixL = int(np.ceil(-Fs*0.5*TB))    	# Left index for time axis
    ixR = int(np.ceil(Fs*(N-0.5)*TB))     # Right index for time axis 
    tt = np.arange(ixL,ixR)/float(Fs)     # Time axis for s(t)
    #print(tt)
    # ***** Conversion from DT a_n to CT a_s(t) *****
    ast = np.zeros(len(tt))    	# Initialize a_s(t)
    ix = array(around(Fs*arange(0,N)*TB),int)    # Symbol center indexes
    ast[ix-int(ixL)] = Fs*an    # delta_n -> delta(t) conversion
    
    # ***** Set up PAM pulse p(t) *****
    ptype = ptype.lower()    # Convert ptype to lowercase
    
    # Set left/right limits for p(t)
    if (ptype=='rect' or ptype=='man'):
        kL = -0.5; kR = -kL
    else:
        kL = -1.0; kR = -kL
    
    # Default left/right limits
    ixpL = np.ceil(Fs*kL*TB)    	# Left index for p(t) time axis
    ixpR = np.ceil(Fs*kR*TB)    	# Right index for p(t) time axis
    ttp =  np.arange(ixpL,ixpR)/float(Fs)     # Time axis for p(t)
    pt = zeros(len(ttp))    	# Initialize pulse p(t)
    if (ptype=='rect'):    	# Rectangular p(t)
        ix = where(logical_and(ttp>=kL*TB, ttp<kR*TB))[0]
        pt[ix] = ones(len(ix))
    elif (ptype == 'tri'):
        pt = array([(1+i*1/TB) if i*1/TB+1 < 1.0 else (1-i*1/TB) for i in list(ttp)])
    elif (ptype=='sinc'):
        k=pparms[0]
        kL = kL*k; kR = -kL
        ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
        ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
        ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
        pt = zeros(len(ttp))
        beta=pparms[1]
        pt = array([sin(pi*t/TB)/(pi*t/TB) if t!= 0  else  1.0 for t in list(ttp)])
        pt=pt*kaiser(len(pt), beta)
    elif (ptype=='man'):
        ix = where(logical_and(ttp>=kL*TB, ttp<=kR*TB))[0]
        pt[ix] = ones(len(ix))
        pt[ix] = concatenate([-1*ones(int(len(ix)/2)),ones(len(pt)-int(len(pt)/2))])
        #for i in ix:
            #print(pt[i])
    
    elif (ptype=='rcf'):
        k=pparms[0]
        alpha=pparms[1]
        kL = kL*k; kR = -kL
        ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
        ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
        ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
        pt = np.pi/4*ones(len(ttp))
        ix0= where(ttp==0)[0]
        pt[ix0] = array([1.0])
        ix0 = where(ttp==TB/(2*alpha))[0]
        pt[ix0] = array([0])
        ix0 = where(ttp==-TB/(2*alpha))[0]
        pt[ix0] = array([0])
        ix_rest = where(logical_and(logical_and(ttp!=0, ttp != TB/(2*alpha)), ttp != -TB/(2*alpha) ))[0]
        for t in ttp:
            num = (sin(pi*ttp/TB)/(pi*ttp/TB))*(cos(pi*alpha*ttp/TB))
            denum =(1-np.power(2*alpha*ttp/TB))
        
    
    elif (ptype=='pr1'):
        k=pparms[0]
        kL = kL*k; kR = -kL
        ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
        ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
        ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
        pt = zeros(len(ttp))
        beta=pparms[1]
        pt = array([sin(pi*t/TB)/(pi*t*(1-t/TB)/TB) if t!= 0 and t!=TB  else  1.0 for t in list(ttp)])
        pt=pt*kaiser(len(pt), beta)  
    elif (ptype=='pr1'):
        k=pparms[0]
        kL = kL*k; kR = -kL
        ixpL = ceil(Fs*kL*TB)		# Left index for p(t) time axis
        ixpR = ceil(Fs*kR*TB)		# Right index for p(t) time axis
        ttp = arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
        pt = zeros(len(ttp))
        beta=pparms[1]
        pt = array([sin(pi*t/TB)/(pi*t*(1-t/TB)/TB) if t!= 0 and t!=TB  else  1.0 for t in list(ttp)])
        pt=pt*kaiser(len(pt), beta)   
    else:
        print("ptype '%s' is not recognized" % ptype)
    
    # ***** Filter with h(t) = p(t) *****
    
    st = convolve(ast,pt)/float(Fs)     # s(t) = a_s(t)*p(t)
    st = st[int(-ixpL):int(ixR-ixL-ixpL)]     	# Trim after convolution
    return comsig.sigWave(st,Fs,0)

