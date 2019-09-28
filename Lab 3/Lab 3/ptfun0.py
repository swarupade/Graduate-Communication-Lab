# File ptfun0.py
# Functions for  gnuradio-companion PAM p(t) generation
import numpy as np

def pampt(sps, ptype, pparms=[]):
    """
    PAM pulse p(t) = p(n*TB/sps), at t=n*Ts, generation
    with 'samples per symbol' sps = TB/Ts = Fs/FB in terms
    of sampling rate Fs = 1/Ts and baud rate FB = 1/TB.
    >>>>> pt = pampt(sps, ptype, pparms) <<<<<
    where  sps:    samples per symbol (Fs/FB)
           ptype:  pulse type ('rect', 'sinc', 'tri')
           pparms  not used for 'rect', 'tri'
           pparms = [k, beta] for 'sinc'
           k:      'tail' truncation parameter for 'sinc',
                   truncates p(t) to -k*sps <= n < k*sps
           beta:   Kaiser window parameter for 'sinc'
           pt:     pulse p(t) at t=n*Ts=n*TB/sps 
    """
    if ptype.lower() == 'rect':
        nn = np.arange(sps)
        pt = np.ones(nn.size)
    else:
        pt = np.ones(1)
    return pt
    
