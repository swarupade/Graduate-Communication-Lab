import numpy as np 
import comsig 
import PAM11 
import showfun1 as showfun

Fs = 44100 # Sampling rate 
FB = 200 # Baud rate 
NTd = 50 # Number of traces displayed 
N = NTd+10 # Number of data symbols 
L = 4 # Number of data levels 
dly = 0.5 # Trigger delay TB/2 
dn = np.floor(L*np.random.rand(N)) # Unipolar L-valued random data 
an = 2*dn - (L-1) # Polar L-valued DT sequence 
sig_an = comsig.sigSequ(an, FB, 0) # DT sequence class 
tt, sig_st = PAM11.pam11(sig_an, FB, Fs, "sinc", [20, 0]) # PAM signal, ’sinc’ p(t) 
showfun.showeye(sig_st, Fs, FB, NTd, [dly, 3, -1.5*L, 1.5*L])
