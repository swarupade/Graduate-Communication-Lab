import numpy as np 
import comsig 
import showfun1  

Fs = 44100 # Sampling rate 
tlen = 1 # Duration in seconds 
tp = 0.1  # Rectangular pulse width
no_ones = int(tp*Fs)  # No of ones in the sequence
no_zeros = int(tlen*Fs - no_ones)    #no_of_zeros
xt = np.array(np.hstack((np.zeros(int(no_zeros/2)), np.ones(no_ones),np.zeros(int(tlen*Fs)-no_ones-int(no_zeros/2)))),np.int8)
sig_xt = comsig.sigWave(xt, Fs, 0) # Combined sinusoidal signal 
N = sig_xt.Nsamp
showfun1.showpsd(sig_xt, Fs,[-40, 41,-60], N) #Plot S_x(f)
