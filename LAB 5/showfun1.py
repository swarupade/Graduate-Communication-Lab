# File: showfun.py
# "show" functions like showft, showpsd, etc
from pylab import *
import numpy as np 
import matplotlib.pyplot as plt
def showft(tt, sig_xt, Fs, ff_lim):
    """
    Plot (DFT/FFT approximation to) Fourier transform of x(t)
    Displays magnitude |X(f)| either linear and absolute or
    normalized (wrt to maximum value) in dB. Phase of X(f) is
        shown in degrees.
        >>>>> showft(tt, xt, ff_lim) <<<<<
        where tt:time axis (increments Ts=1/Fs) for x(t)
            xt:sampled CT signal x(t))
            ff_lim = [f1,f2,llim]
            f1:lower frequency limit for display
            f2:upper frequency limit for display
            llim = 0: display |X(f)| linear and absolute
            llim > 0: same as llim = 0 but phase is masked
                      (set to zero) for |X(f)| < llim
            llim < 0: display 20*log_{10}(|X(f)|/max(|X(f)|))
                      in dB with lower display limit llim dB,
                      phase is masked (set to zero) for f
                      with magnitude (dB, normalized) < llim
        """
    # ***** Prepare x(t), swap pos/neg parts of time axis *****
    n = len(tt)
    ixp = np.where(tt>=0)[0]			# Indexes for t>=0
    ixn = np.where(tt<0)[0]			# Indexes for t<0
    tlen = tt[-1]	
    xt = sig_xt.signal()
    xt = np.hstack((xt[ixp],xt[ixn])) 		# Swap pos/neg time axis parts
    llim = ff_lim[2]
    f11, f12 = ff_lim[0], ff_lim[1]
    
    
    # ***** Compute X(f), make frequency axis *****
    Xf = np.fft.fft(xt)/float(Fs)			# DFT/FFT of x(t),
    
    # scaled for X(f) approximation
    ff = Fs*np.arange(n)/float(n) 		# Frequency axis
    if ff_lim[0] < 0:
    	ixp = np.where(ff < Fs/2)[0]
    	ixn = np.where(ff >= Fs/2)[0]
    	ff = np.hstack((ff[ixn]-Fs,ff[ixp]))
    	Xf = np.hstack((Xf[ixn],Xf[ixp]))

    # ***** Compute |X(f)|, arg[X(f)] *****
    absXf = abs(Xf)				# Magnitude |X(f)|
    argXf = np.angle(Xf)			# Phase arg[X(f)]
    
    if llim < 0:
    	absXf = 20*np.log10((absXf)/max(absXf))
    	ix = np.where(absXf < ff_lim[2])[0]
    	argXf[ix] = np.zeros(len(ix))
    	absXf[ix] = ff_lim[2]*np.ones(len(ix))
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
    strgt = 'FT Approximation, $F_s=$' + str(Fs) + ' Hz'
    strgt = strgt + ', N=' + str(n)
    strgt = strgt + ', $\Delta_f$={0:3.2f}'.format(Fs/float(n)) + ' Hz'
    af11.set_title(strgt)
    af12 = f1.add_subplot(212)
    af12.plot(ff,(180/pi)*argXf) 		# Plot phase in degrees
    af12.grid()
    af12.set_ylabel('arg[X(f)] [deg]')
    af12.set_xlabel('f [Hz]')
    plt.show()

def showeye(xt, Fs, FB, NTd=50, dispparms=[]):
    """
    Display eye diagram of digital PAM signal r(t)
    >>>>> showeye(rt, Fs, FB, NTd, dispparms) <<<<<
    where rt: received PAM signal r(t)=sum_n a_n*q(t-nTB)
    Fs: sampling rate for r(t)
    FB: Baud rate of DT sequence a_n, TB = 1/FB
    NTd: Number of traces displayed
    dispparms = [delay, width, ylim1, ylim2]
    delay: trigger delay (in TB units, e.g., 0.5)
    width: display width (in TB units, e.g., 3)
    ylim1: lower display limit, vertical axis
    ylim2: upper display limit, vertical axis
    """
    rt = xt.signal()
    N = xt.Nsamp
    t0 = dispparms[0]/float(FB)         # Delay in sec
    tw = dispparms[1]/float(FB)         # Display width in sec
    dws = int(np.floor(Fs*tw))                  # Display width in samples
    tteye = np.arange(dws)/float(Fs)       # Time axis for eye
    trix = np.array(np.around(Fs*(t0+np.arange(NTd)/float(FB))),int)
    ix = np.where(np.logical_and(trix>=0, trix<= N-dws))[0]
    trix = trix[ix]                     # Trigger indexes within r(t)
    TM = rt[trix[0]:trix[0]+dws]        # First trace
    #TM = np.vstack((TM, rt[trix[1]:trix[1]+dws]))
    for ind in np.arange(1, NTd):
        TM = np.vstack((TM, rt[trix[ind]:trix[ind]+dws]))    # Second trace
    #plt.figure()
    plt.plot(FB*tteye, TM.T, '-b')
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.title("Eye Diagram for r(t) with FB = %s Baud, \n t0 = %s , # Traces=%s" %(FB,dispparms[0],NTd))
    plt.xlabel('t/tb')
    plt.ylabel('r(t)')
    plt.ylim([dispparms[2],dispparms[3]])
    plt.grid()
    plt.show()

def showpsd(sig_xt, Fs, ff_lim, N):
    """
    Plot (DFT/FFT approximation to) power spectral density (PSD) of x(t).
    Displays S_x(f) either linear and absolute or normalized in dB.
    >>>>> showpsd(sig_xt, Fs, ff_lim, N) <<<<<
    where sig_xt: waveform from class sigWave 
        sig_xt.signal(): sampled CT signal x(t) 
        sig_xt.get_Fs(): sampling rate of x(t) 
        Fs:sampling rate of x(t)
        ff_lim = [f1,f2,llim]
        f1: lower frequency limit for display
        f2:upper frequency limit for display
        llim = 0: display S_x(f) linear and absolute
        llim < 0: display 10*log_{10}(S_x(f))/max(S_x(f))
        in dB with lower display limit llim dB
        N:blocklength
    """
    # ***** Determine number of blocks, prepare x(t) *****

    rt = sig_xt.signal()                               # Get x(t) 
    Fs = Fs                               # Sampling rate of x(t) 
    #xt = diff(hstack((0,rt)))*Fs
    xt = rt*rt*rt*rt
    #xt = rt
    N = int(min(N, len(xt)))                           # N <= length(xt) needed
    NN = int(np.floor(len(xt)/float(N)))                  # Number of blocks of length N
    xt = xt[0:N*NN]                                    # Truncate x(t) to NN blocks
    xNN = np.reshape(xt,(NN,N))                           # NN row vectors of length N
    llim = ff_lim[2]                                   # to display Sxf in db or linear
   
    # ***** Compute DFTs/FFTs, average over NN blocks *****
    Sxf = np.power(abs(fft(xNN)),2.0)                  # NN FFTs, mag squared
    if NN > 1:
        Sxf = sum(Sxf, axis=0)/float(NN)
    
    Sxf = Sxf/float(N*Fs)                               # Correction factor DFT -> PSD
    Sxf = np.reshape(Sxf,size(Sxf))
    ff = Fs*np.array(np.arange(N),int64)/float(N)             # Frequency axis
    
    if ff_lim[0] < 0:                                   # Negative f1 case
        ixp = np.where(ff<0.5*Fs)[0]                       # Indexes of pos frequencies
        ixn = np.where(ff>=0.5*Fs)[0]                      # Indexes of neg frequencies
        ff = np.hstack((ff[ixn]-Fs,ff[ixp]))               # New freq axis
        Sxf = np.hstack((Sxf[ixn],Sxf[ixp]))               # Corresponding S_x(f)
    df = Fs/float(N)                                    # Delta_f, freq sample spacing
    Px = (cumsum(Sxf)*Fs/N)[-1]                         # Calculating total power of the signal
    
    # ***** Determine maximum, trim to ff_lim *****
    maxSxf = max(Sxf)                                   # Maximum of S_x(f)
    
    ixf = np.where(logical_and(ff>=ff_lim[0], ff<ff_lim[1]))[0]
    ff = ff[ixf]                                        # Trim to ff_lim specs
    Sxf = Sxf[ixf]
    strgy = 'PSD: $S_x(f)$'                              # ylabel string
    Px_prime = (cumsum(Sxf)*Fs/N)[-1]                   # Calculating the power of the signal in 
                                                        # frequency range of [ff_lim[0], ff_lim[1]] 
    Px_prime_perc = (Px_prime/Px) * 100.0               # Changing the power content in percentage
    
    # ***** Changing Sxf into dB when llim is less than zero *****
    if llim < 0:
            ixz = where(Sxf==0)
            Sxf[ixz] = 1e-10*ones(len(ixz))
            Sxf  = 10*log10((Sxf/maxSxf))
            ix = where(Sxf <= llim)
            Sxf[ix] = llim*ones(len(ix))
    

    
    # ***** Plot PSD *****
    f1 = plt.figure()
    af1 = f1.add_subplot(111)
    
    if llim < 0:
        strgt = '$P_x=${:0.2f}, $P_x(f1,f2)=${:0.2f}%, $F_s=${:d} Hz\n'.format(Px,Px_prime_perc,Fs)
        strgt1 = 'Sx(f) in dB'
        af1.set_ylim(llim,0)
    else:
        strgt = '$P_x=${:0.2f}, $P_x(f1,f2)=${:0.5f}, $F_s=${:d} Hz\n'.format(Px,Px_prime,Fs)
        strgt1 = 'Sx(f)'
    
    
    strgt = strgt + ', $\\Delta_f=${:.3g} Hz'.format(df)
    strgt = strgt + ', $NN=${:d}, $N=${:d}'.format(NN, N)
    
    af1.plot(ff, Sxf, '-b')
    af1.grid()
    af1.set_xlabel('f [Hz]')
    af1.set_ylabel(strgt1)
    af1.set_title(strgt)
    plt.show()