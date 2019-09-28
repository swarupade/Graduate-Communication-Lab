# File: showfun.py 
# "show" functions like showft, showpsd, etc 
from pylab import *
import numpy as np 
import matplotlib.pyplot as plt
#import ftpam01 as ft
def showft(tt, sig_xt, Fs, ff_lim): 
   
    """ Plot (DFT/FFT approximation to) Fourier transform of waveform x(t). 
    Displays magnitude |X(f)| either linear and absolute or normalized (wrt to maximum value) in dB. 
    Phase of X(f) is shown in degrees. 
    >>>>> showft(sig_xt, ff_lim) <<<<< 
    where sig_xt: waveform from class sigWave 
    sig_xt.timeAxis(): time axis for x(t) 
    sig_xt.signal(): sampled CT signal x(t) 
    ff_lim = [f1, f2, llim] 
    f1: lower frequency limit for display 
    f2: upper frequency limit for display 
    llim = 0: display |X(f)| linear and absolute 
    llim > 0: same as llim = 0 but 
              phase is masked (set to zero) for |X(f)| < llim 
    llim < 0: display 20*log_{10}(|X(f)|/max(|X(f)|)) 
              in dB with lower display limit llim dB, 
              phase is masked (set to zero) for f with magnitude (dB, normalized) 
              < llim 
    """
    N = len(tt) # Blocklength of DFT/FFT 
    ixp = np.where(tt>=0)[0] # Indexes for t>=0 
    ixn = np.where(tt<0)[0] # Indexes for t<0 
    tlen = tt[-1]
    xt = sig_xt.signal() # Get x(t) 
    xt = np.hstack((xt[ixp],xt[ixn]))
    llim = ff_lim[2]
    f11, f12 = ff_lim[0], ff_lim[1]

    # ***** Compute X(f), make frequency axis ***** 

    Xf = np.fft.fft(xt)/float(Fs) # DFT/FFT of x(t), # scaled for X(f) approximation
    ff = Fs*np.arange(N)/float(N) # Frequency axis 
    if ff_lim[0] < 0:
        ixp = np.where(ff < Fs/2)[0]
        ixn = np.where(ff >= Fs/2)[0]
        ff = np.hstack((ff[ixn]-Fs,ff[ixp]))
        Xf = np.hstack((Xf[ixn],Xf[ixp]))

    # ***** Compute |X(f)|, arg[X(f)] ***** 
    absXf = np.abs(Xf)				# Magnitude |X(f)|
    argXf = np.angle(Xf)			# Phase arg[X(f)]
     # ***** Mirror |X(f)| about 0 (if ff_lim[0]<0) *****
    if f11<0:
        absXf = concatenate([absXf[::-1],absXf])
        neg = [-1*i for i in argXf]
        argXf = concatenate([neg[::-1],argXf])
        neg = [-1*i for i in ff]
        ff = concatenate([neg[::-1],ff])
    # ***** Floor values of argXf for points where absXf<llim *****
    if llim>0:
        for i in range(0,len(absXf)):
            if absXf[i] < llim:
                argXf[i] = 0
    # ***** Convert absXt to dB and floor argXf for points where absXf<llim(dB) *****
    if llim<0:
        mag=10**(llim/20)
        absXfmax=amax(absXf)
        for i in range(0,len(absXf)):
            if absXf[i]>mag:
                absXf[i] = 20*math.log10(absXf[i]/absXfmax)
            else:
                absXf[i]=llim
                argXf[i]=0
    
   # ***** Plot magnitude/phase *****
    f1 = plt.figure()
    af11 = f1.add_subplot(211)
    af11.plot(ff,absXf)
    af11.grid()
    if llim >= 0:
        af11.set_ylabel('|X(f)|')
    else:
	    af11.set_ylabel('|X(f)| in dB')
    strgt = 'FT Approximation, $F_s=10000'  + ' Hz'
    strgt = strgt + ', N=' + str(N)
    strgt = strgt + ', $\Delta_f$={0:3.2f}'.format(10000/float(N)) + ' Hz'
    af11.set_title(strgt)
    xlim([ff_lim[0],ff_lim[1]])
    af12 = f1.add_subplot(212)
    af12.plot(ff,(180/np.pi)*argXf) 		# Plot phase in degrees
    af12.grid()
    af12.set_ylabel('arg[X(f)] [deg]')
    af12.set_xlabel('f [Hz]')
    xlim([ff_lim[0],ff_lim[1]])
    show()
    
    

