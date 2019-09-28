import numpy as np 
import comsig 
import showfun1  
import PAM11

tlen  = 5
Fb = 100
Fs = 44100
no_bits = int(tlen*Fb)
an = 2*(np.around(2*np.random.rand(no_bits)) %2 ) - 1
sig_an = comsig.sigWave(an,Fs,0)
tt, st = PAM11.pam11(sig_an, Fb, Fs, 'sinc',[20,0])
N = Fs
showfun1.showpsd(st, Fs,[-200,200, -60], N)