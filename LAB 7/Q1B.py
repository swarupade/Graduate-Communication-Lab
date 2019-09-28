import filtfun
import amfun
import showfun1
import numpy as np
import comsig

Fs = 44100 # Sampling rate 
tlen = 1.0 # Duration 
f0, f1 = 3000, 5000 # Message frequencies 
tt = np.arange(np.round(tlen*Fs))/float(Fs) # Time axis 
mt = np.cos(2*np.pi*f0*tt) + np.cos(2*np.pi*f1*tt) # Message signal 
sig_mt = comsig.sigWave(mt, Fs, 0) # Waveform from class sigWave
xtype = 'tc'
fcparms =[9000, 0,0.05]
fmparms =[4000,10,0.7]
sig_xt = amfun.amxmtr(tt,sig_mt,xtype,fcparms,fmparms)
showfun1.showpsd(sig_xt, Fs, [0,18000, 0],sig_xt.Nsamp)