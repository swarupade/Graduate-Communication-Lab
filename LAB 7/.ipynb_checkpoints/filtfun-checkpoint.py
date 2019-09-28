# File: filtfun.py
# Module for filter functions

from pylab import *
from scipy.signal import butter, lfilter
import comsig


def trapfilt(sig_xt, fparms, k, alfa, ftype):
    """
    Delay compensated FIR LPF/BPF filter with trapezoidal
    frequency response.
    >>>>> sig_yt, n = trapfilt(sig_xt, fparms, k, alfa) <<<<<
    where sig_yt: waveform from class sigWave
          sig_yt.signal(): filter output y(t), samp rate Fs
          n: filter order
          sig_xt: waveform from class sigWave
          sig_xt.signal(): filter input x(t), samp rate Fs
          sig_xt.get_Fs(): sampling rate for x(t), y(t)
          fparms = fL for LPF
          fL: LPF cutoff frequency (-6 dB) in Hz
          fparms = [fBW, fc] for BPF
          fBW: BPF -6dB bandwidth in Hz
          fc: BPF center frequency in Hz
          k: h(t) is truncated to
          |t| <= k/(2*fL) for LPF
          |t| <= k/fBW for BPF
          alfa: frequency rolloff parameter, linear
                rolloff over range
          (1-alfa)fL <= |f| <= (1+alfa)fL for LPF
          (1-alfa)fBW/2 <= |f| <= (1+alfa)fBW/2 for BPF
"""
    xt = sig_xt.signal() # Input signal
    Fs = sig_xt.get_Fs() # Sampling rate
    if ftype == 'LPF':
        fL = fparms[0]
        fc = fparms[1]
    elif ftype == 'BPF':
        fBW = fparms[0]
        fc  = fparms[1]
        fL  = fBW/2
   
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
    if fc != 0:
        ht = 2*ht*cos(2*pi*fc*tth)
   
    yt = lfilter(ht, 1, hstack((xt, zeros(ixk))))/float(Fs) # Compute filter output y(t)
    yt = yt[ixk:] # Filter delay compensation
    return comsig.sigWave(yt*Fs,Fs,-0.5),n # Return y(t) and filter order

def trapfilt_cc(sig_xt, fparms, k, alfa, ftype):
    """
      Delay compensated FIR LPF/BPF filter with trapezoidal
      frequency response, complex-valued input/output and
      complex-valued filter coefficients.
    >>>>> sig_yt, n = trapfilt_cc(sig_xt, fparms, k, alfa) <<<<<
          where sig_yt: waveform from class sigWave
          sig_yt.signal(): complex filter output y(t), samp rate Fs
          n: filter order
          sig_xt: waveform from class sigWave
          sig_xt.signal(): complex filter input x(t), samp rate Fs
          sig_xt.get_Fs(): sampling rate for x(t), y(t)
          fparms = fL for LPF
          fL: LPF cutoff frequency (-6 dB) in Hz
          fparms = [fBW, fBc] for BPF
          fBW: BPF -6dB bandwidth in Hz
          fBc: BPF center frequency (pos/neg) in Hz
          k: h(t) is truncated to
          |t| <= k/(2*fL) for LPF
          |t| <= k/fBW for BPF
          alfa: frequency rolloff parameter, linear
               rolloff over range
         (1-alfa)*fL <= |f| <= (1+alfa)*fL for LPF
         (1-alfa)*fBW/2 <= |f| <= (1+alfa)*fBW/2 for BPF
     """
    xt = sig_xt.signal() # Input signal
    Fs = sig_xt.get_Fs() # Sampling rate
    if ftype == 'LPF':
        fL = fparms[0]
        fc = fparms[1]
    elif ftype == 'BPF':
        fBW = fparms[0]
        fc  = fparms[1]
        fL  = fBW/2
    ixk = round(Fs*k/float(2*fL))                                # Tail cutoff index
    tt = arange(-ixk,ixk+1)/float(Fs)                            # Time axis for h(t)
    n = len(tt)-1                                                # Filter order
    # ***** Generate impulse response ht here *****
    ht = zeros(len(tt))                                          # Initializing ht
    ix = where(tt != 0)[0]
    if alfa != 0:
        ht[ix] = ((sin(2*pi*fL*tt[ix]))/(pi*tt[ix]))*((sin(2*pi*alfa*fL*tt[ix]))/(2*pi*alfa*fL*tt[ix]))
    else:
        ht[ix] = (sin(2*pi*fL*tt[ix]))/(pi*tt[ix])
    ix0 = where(tt == 0)[0]
    ht[ix0] = 2*fL
    if fc != 0:
        ht = 2*ht*exp(1j*2*pi*fc*tt)
    yt = lfilter(ht, 1, hstack((xt, zeros(ixk))))/float(Fs)      # Compute filter output y(t)
    yt = yt[ixk:]                                                # Filter delay compensation
    return comsig.sigWave(yt, Fs, -0.5), n                                                 # Return y(t) and filter order


