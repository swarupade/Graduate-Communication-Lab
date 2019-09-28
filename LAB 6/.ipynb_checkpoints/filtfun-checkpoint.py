# File: filtfun.py
# Module for filter functions

from pylab import *
from scipy.signal import butter, lfilter
import comsig


def trapfilt(sig_xt, fL, k, alfa):
    """
    Delay compensated FIR LPF/BPF filter with trapezoidal frequency response.
    >>>>> sig_yt, n = trapfilt(sig_xt, fL, k, alfa) <<<<<
    where
        sig_yt: waveform from class sigWave
            sig_yt.signal(): filter output y(t), samp rate Fs
            n: filter order
        sig_xt: waveform from class sigWave
            sig_xt.signal(): filter input x(t), samp rate Fs
            sig_xt.get_Fs(): sampling rate for x(t), y(t)
        fL: LPF cutoff frequency (-6 dB) in Hz
        k: h(t) is truncated to |t| <= k/(2*fL)
        alfa: frequency rolloff parameter, linear rolloff over range (1-alfa)fL <= |f| <= (1+alfa)fL
    """
    xt = sig_xt.signal() # Input signal
    Fs = sig_xt.get_Fs() # Sampling rate
    ixk = round(Fs*k/float(2*fL)) # Tail cutoff index
    tth = arange(-ixk,ixk+1)/float(Fs) # Time axis for h(t)
    n = len(tth)-1 # Filter order
    ht_num = (sin(2*pi*fL*tth)*sin(2*pi*alfa*fL*tth))
    ht_den = (pi*tth*2*pi*alfa*fL*tth)
    nans = where(ht_den==0)
    for i in nans:
        ht_den[i]=1 # ht[i-1]
        ht_num[i]=2*fL # ht[i-1]
    ht = ht_num/ht_den
    titlestr = 'Trapezoidal LPF, $h_L(t)$ Truncated to $|t|<k/(2f_L)$, $f_L$ = '+str(fL)+' Hz, k = '+str(k)+', alpha = '+str(alfa)
    yt = lfilter(ht, 1, hstack((xt, zeros(ixk))))/float(Fs) # Compute filter output y(t)
    yt = yt[ixk:] # Filter delay compensation
    return ecen.sigWave(yt*Fs, Fs, sig_xt.get_t0()), n # Return y(t) and filter order