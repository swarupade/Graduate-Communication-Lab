import amfun
import comsig
import showfun1
import wavfun
import matplotlib.pyplot as plt
import numpy as np

Fs = 44100

FB, rt1 = wavfun.wavread('myam701.wav')
FB, rt2 = wavfun.wavread('myam702.wav')
FB, rt3 = wavfun.wavread('myam703.wav')
rt3_sig = comsig.sigWave(rt3, Fs, 0)
rt1_sig = comsig.sigWave(rt1, Fs, 0)
rt2_sig = comsig.sigWave(rt2, Fs, 0)
fcparms1 =[8000, (-90*np.pi)/180]
fcparms2 =[8000, 0]
fmparms  =[4000,10,0.05]
sig_mthat1 = amfun.amrcvr(rt3_sig.timeAxis(),rt3_sig,'coh',fcparms1, fmparms,[],'LPF')
sig_mthat2 = amfun.amrcvr(rt3_sig.timeAxis(),rt3_sig,'coh',fcparms2, fmparms,[],'LPF')
wavfun.wavwrite('recovered703speech.wav', Fs, 0.999*sig_mthat1.signal()/float(max(abs(sig_mthat1.signal()))))
wavfun.wavwrite('recovered702music.wav', Fs, 0.999*sig_mthat2.signal()/float(max(abs(sig_mthat2.signal()))))


