import wavfun
import comsig
import showfun1
import filtfun
from pylab import *

xt, framerate = wavfun.wavread("pamsig601.wav")
sig_xt = comsig.sigWave(xt, framerate)
#showfun1.showpsd(sig_xt, framerate, [-1000,1000,0],framerate)

fL=1000
k=10
alfa=0.2
[trap_xt,n]=filtfun.trapfilt(sig_xt,fL,k,alfa)
trap_xt_sqrd = trap_xt.copy()

trap_xt_sqrd.sig = trap_xt_sqrd.sig**2
#showfun1.showpsd(trap_xt_sqrd,framerate,[-1000, 1000, 0],framerate) #Plot S_x(f)
#showfun1.showpsd(trap_xt_sqrd,framerate,[10, 1000, 0],framerate) #Plot S_x(f)
showfun1.showpsd(trap_xt_sqrd,framerate,[810, 820, 0],framerate) #Plot S_x(f)