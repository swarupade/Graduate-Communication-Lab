from pylab import *
import comsig
import showfun1
import filtfun

Fs=16000
t0=-0.5
t=1
xt = concatenate([zeros(int(Fs/2)),[1],zeros(Fs-int(Fs/2))])

sig_xt= comsig.sigWave(xt, Fs, t0)
fL=1000
k=20
alfa=0

[trapfilt_xt, FS]=filtfun.trapfilt(sig_xt,fL,k,alfa)
tt = trapfilt_xt.timeAxis()
showfun1.showft(tt,trapfilt_xt,Fs,[-3000,3000,0])
showfun1.showft(tt,trapfilt_xt,Fs,[-3000,3000,-60])