from pylab import *
import comsig
import matplotlib.pyplot as plt
import ascfun
import PAM12


Fb = 100 # Baud rate (bits/sec)
Fs = 44100 # Sampling frequency (samples/sec)

string="Test"
dn = ascfun.asc2bin(string,8)
dn = multiply(dn,2)-1

sig_pt = comsig.sigSequ(dn,Fb)
tt, pam_pt = PAM12.pam12(sig_pt,Fb, Fs,'rrcf',[5,0.4])
plt.plot(tt,pam_pt.signal())