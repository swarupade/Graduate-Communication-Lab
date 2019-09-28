# File: filtspecs.py
# Module for filter specfifications
# FIR: Determine filter taps
# IIR: Determine numerator (b) and denominator (a)
# polynomial coefficients

from pylab import *
from scipy.signal import butter

def trapfilt_taps(N, phiL, alfa):
    """
    Returns taps for order N FIR LPF with trapezoidal frequency
    response, normalized cutoff frequency phiL = fL/Fs, and rolloff
    parameter alfa.
    >>>>> hLn = trapfilt_taps(N, phiL, alfa) <<<<<
    where
        N: filter order
        phiL: normalized cutoff frequency (-6 dB)
        alfa: frequency rolloff parameter, linear rolloff over range (1-alfa)phiL <= |f| <= (1+alfa)phiL
    """
    tth = arange(-N/2.0,N-(N/2.0)) # Time axis for h(t)
    hLn_num = (sin(2*pi*phiL*tth)*sin(2*pi*alfa*phiL*tth))
    hLn_den = (pi*tth*2*pi*alfa*phiL*tth)
    nans = where(hLn_den==0)
    for i in nans:
        hLn_den[i]=1 # ht[i-1]
        hLn_num[i]=2*phiL # ht[i-1]
    hLn = hLn_num/hLn_den
    return hLn