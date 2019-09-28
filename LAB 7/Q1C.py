import filtfun
import amfun
import showfun1
import numpy as np
import comsig
import wavfun

Fs = 44100 # Sampling rate 
#speech and music signal for modulating
FB1, mt1 = wavfun.wavread('speech701.wav')
FB2, mt2 = wavfun.wavread('music701.wav')
sig_mt1 = comsig.sigWave(mt1, Fs, 0) # Waveform from class sigWave
sig_mt2 = comsig.sigWave(mt2, Fs, 0)
sig_mt2ad = comsig.sigWave(0.7*sig_mt2.signal(), Fs)
# generating  AM-DSB-SC signals

xtype = 'sc'
fcparms1 =[8000,(-90*np.pi)/180]
fcparms2 =[8000, 0]
fmparms  =[4000,10,0.05]
sig_xt1 = amfun.amxmtr(sig_mt1.timeAxis(),sig_mt1,xtype,fcparms1,fmparms)

sig_xt2 = amfun.amxmtr(sig_mt2ad.timeAxis(),sig_mt2ad,xtype,fcparms2,fmparms)
#ShowPSD for the first two

showfun1.showpsd(sig_xt1, Fs, [-18000,18000, 0],sig_xt1.Nsamp)
showfun1.showpsd(sig_xt2, Fs, [-18000,18000, 0],sig_xt2.Nsamp)

#The third signal
sig_xt3 = comsig.sigWave((sig_xt1.signal()+sig_xt2.signal())/np.sqrt(2),Fs)
showfun1.showpsd(sig_xt3, Fs, [-18000,18000, 0],sig_xt3.Nsamp)

#Writing the signals to wav files

wavfun.wavwrite('myam701.wav',Fs,(0.999*sig_xt1.signal()/float(max(abs(sig_xt1.signal())))))
wavfun.wavwrite('myam702.wav',Fs,(0.999*sig_xt2.signal()/float(max(abs(sig_xt2.signal())))))
wavfun.wavwrite('myam703.wav',Fs,(0.999*sig_xt3.signal()/float(max(abs(sig_xt3.signal())))))

