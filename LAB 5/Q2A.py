import comsig 
import numpy as np
import showfun1 
import PAM11
import wavfun

Fs = 44100 # Sampling rate
nFl = 1000
tlen = 5 # Duration in seconds 
tt = np.arange(np.round(tlen*Fs))/float(Fs) # Time axis 
nn = np.random.randn(np.round(tlen * 2 * nFl)) # Gaussian noise n(t) 
sig_nn = comsig.sigSequ(nn, 2*nFl, 0)
N = Fs
tt, sig_nt = PAM11.pam11(sig_nn, 2*nFl , Fs,'rcf',[20, 0.2]) 
showfun1.showpsd(sig_nt, Fs, [-10000,10000,-60],N)
wavfun.wavwrite('bandwidth1000.wav', Fs, sig_nt.signal())

