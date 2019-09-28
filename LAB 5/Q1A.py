import numpy as np 
import comsig 
import showfun1  
Fs = 44100 # Sampling rate 
f1 = 700 # Test frequency 1 
f2 = 720 # Test frequency 1 
tlen = 2 # Duration in seconds 
tt = np.arange(round(tlen*Fs))/float(Fs) # Time axis 
x1t = np.sin(2*np.pi*f1*tt) # Sine with freq f1 
x2t = 0.01*np.cos(2*np.pi*f2*tt) # Attenuated cosine with freq f2 
sig_xt = comsig.sigWave(x1t+x2t, Fs, 0) # Combined sinusoidal signal 
N = sig_xt.Nsamp
showfun1.showpsd(sig_xt, Fs,[-1000, 1000, -100], N) #Plot S_x(f)
