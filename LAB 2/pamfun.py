# File: pamfun.py 
# Functions for pulse amplitude modulation (PAM) import numpy as np import comsig
import numpy as np
import comsig
def pam10(sig_an, Fb, Fs, ptype, pparms=[]): 
    """ Pulse amplitude modulation: 
    a_n -> s(t), -TB/2<=t<(N-1/2)*TB, 
    V1.0 for ’rect’, ’sinc’, and ’tri’ pulse types. >>>>> 
    sig_st = pam10(sig_an, Fs, ptype, pparms) <<<<< 
    where sig_an: sequence from class sigSequ 
    sig_an.signal(): N-symbol DT input sequence a_n, 0 <= n < N 
    sig_an.get_FB(): Baud rate of a_n, TB=1/FB 
    Fs: sampling rate of s(t) 
    ptype: pulse type (’rect’,’sinc’,’tri’) 
    pparms not used for ’rect’,’tri’ pparms = [k, beta] for ’sinc’ 
    k: "tail" truncation parameter for ’sinc’ (truncates p(t) to -k*TB <= t < k*TB) 
    beta: Kaiser window parameter for ’sinc’ 
    sig_st: waveform from class sigWave 
    sig_st.timeAxis(): time axis for s(t), starts at -TB/2 
    sig_st.signal(): CT output signal s(t), -TB/2<=t<(N-1/2)*TB, with sampling rate Fs 
    """
    N = len(sig_an) # Number of data symbols 
    FB = Fb # Baud rate 
    TB = 1/FB
    n0 = 0# Starting index 
    ixL = int(np.ceil(Fs*(n0-0.5)/float(FB))) # Left index for time axis 
    ixR = int(np.ceil(Fs*(n0+N-0.5)/float(FB))) # Right index for time axis 
    tt = np.arange(ixL,ixR)/float(Fs) # Time axis for s(t) 
    t0 = tt[0] # Start time for s(t) 
    
    # ***** Conversion from DT a_n to CT a_s(t) ***** 
    an = sig_an # Sequence a_n 
    ast = np.zeros(len(tt)) # Initialize a_s(t) 
    ix = np.array(np.around(Fs*(np.arange(0,N)+n0)/float(FB)),int) # Symbol center indexes 
    ast[ix-int(ixL)] = Fs*an # delta_n -> delta(t) onversion 
    # ***** Set up PAM pulse p(t) ***** 
    ptype = ptype.lower() # Convert ptype to lowercase
                          # Set left/right limits for p(t) 
    if (ptype=="rect"): 
        kL = -0.5; kR = -kL 
    else: 
        kL = -1.0; kR = -kL # Default left/right limits 
    ixpL = int(np.ceil(Fs*kL/float(FB))) # Left index for p(t) time axis 
    ixpR = int(np.ceil(Fs*kR/float(FB))) # Right index for p(t) time axis 
    ttp = np.arange(ixpL,ixpR)/float(Fs) # Time axis for p(t) 
    pt = np.zeros(len(ttp)) # Initialize pulse p(t) 
    if (ptype=="rect"): # Rectangular p(t) 
        ix = np.where(np.logical_and(ttp>=kL/float(FB),ttp<kR/float(FB)))[0] 
        pt[ix] = np.ones(len(ix))
    elif (ptype == "tri"):
	    pt = np.array([(1+i*1/TB) if i*1/TB+1 < 1.0 else (1-i*1/TB) for i in list(ttp)]) 
    elif (ptype=="sinc"):
        k=pparms[0]
        kL = -1.0*k
        kR = -kL
        ixpL = int(np.ceil(Fs*kL*TB)) 		# Left index for p(t) time axis
        ixpR = int(np.ceil(Fs*kR*TB))		# Right index for p(t) time axis
        ttp = np.arange(ixpL,ixpR)/float(Fs)	 # Time axis for p(t)
        pt = np.zeros(len(ttp))
        beta=pparms[1]
        pt = np.array([np.sin(np.pi*t/TB)/(np.pi*t/TB) if t!= 0  else  1.0 for t in list(ttp)])
        pt=pt*np.kaiser(len(pt), beta)
    else:
        print("ptype ’%s’ is not recognized" % ptype)
    # ***** Filter with h(t) = p(t) ***** 
    st = np.convolve(ast, pt)/float(Fs)  # s(t) = a_s(t)*p(t) 
    st = st[-ixpL:ixR-ixL-ixpL] # Trim after convolution
    return(tt,st)
    return comsig.sigWave(st, Fs, t0) # Return waveform from sigWave class
