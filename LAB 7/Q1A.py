from pylab import *
import comsig
import showfun1
import filtfun


Fs=44100
t0=-0.5
t=1
xt=concatenate([zeros(int(Fs/2)),[1],zeros(Fs-int(Fs/2))])

sig_xt=comsig.sigWave(xt, Fs, t0)

fparms = [7000,10500]
k=20
alfa= 0.05
ftype = 'BPF'
[trapfilt_xt, n]=filtfun.trapfilt(sig_xt,fparms,k,alfa,ftype)
tt = trapfilt_xt.timeAxis()
showfun1.showft(tt,trapfilt_xt,Fs,[0,15000,0])
#showfun1.showft(tt,trapfilt_xt,Fs,[0,15000,-60])