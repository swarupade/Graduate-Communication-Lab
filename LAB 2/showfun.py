# File: showfun.py 
# "show" functions like showft, showpsd, etc 
import numpy as np 
import matplotlib.pyplot as plt
import ftpam01 as ft
def showft(tt,sig_xt, Fs, ff_lim): 
   
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
    N = len(sig_xt) # Blocklength of DFT/FFT 
    ixp = np.where(tt>=0)[0] # Indexes for t>=0 
    ixn = np.where(tt<0)[0] # Indexes for t<0 
    xt = sig_xt # Get x(t) 
    xt = np.hstack((xt[ixp],xt[ixn]))
    llim = ff_lim[2]
    f11, f12 = ff_lim[0], ff_lim[1]

    # ***** Compute X(f), make frequency axis ***** 

    Xf = np.fft.fft(xt)/float(Fs) # DFT/FFT of x(t), # scaled for X(f) approximation
    ff = Fs*np.array(np.arange(N),np.int64)/float(N) # Frequency axis 
    if ff_lim[0] < 0:
        ixp = np.where(ff < Fs/2)[0]
        ixn = np.where(ff >= Fs/2)[0]
        ff = np.hstack((ff[ixn]-Fs,ff[ixp]))
        Xf = np.hstack((Xf[ixn],Xf[ixp]))

    # ***** Compute |X(f)|, arg[X(f)] ***** 
    absXf = np.abs(Xf)				# Magnitude |X(f)|
    argXf = np.angle(Xf)			# Phase arg[X(f)]
	
    if llim < 0:
        absXf = 20*np.log10((absXf)/max(absXf))
        ix = np.where(absXf < ff_lim[2])[0]
        argXf[ix] = np.zeros(len(ix))
        absXf[ix] = -60*np.ones(len(ix))
    else:
	    ix = np.where(absXf < ff_lim[2])[0]
	    argXf[ix]=np.zeros(len(ix))
        
    ix = np.where(np.logical_and(ff>= ff_lim[0], ff < ff_lim[1]))[0]
    ff = ff[ix]
    absXf=absXf[ix]
    argXf=argXf[ix]
    
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
    af12 = f1.add_subplot(212)
    af12.plot(ff,180/np.pi*argXf) 		# Plot phase in degrees
    af12.grid()
    af12.set_ylabel('arg[X(f)] [deg]')
    af12.set_xlabel('f [Hz]')

