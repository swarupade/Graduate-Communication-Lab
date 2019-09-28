#plotting of the signal

import wavfun
import sinc_ipol
from pylab import *
from numpy import *

#from ast import literal_eval
 
fO = 100
tlen = 5e-2
Fs, rt = wavfun.wavread("sig01.wav")
print(Fs)
tt = arange(0,round(tlen*Fs))/float(Fs)
plot(tt[0:60],rt[0:60]), grid()
title("Original Signal")
show()

r3at=vstack([rt,zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt)),zeros(len(rt))])
r3t = reshape(r3at,r3at.size,order="F")
fs = Fs*120
t3t = arange(0,len(r3t))/float(fs)
#print(tt[0:150])
#print(rt[0:150])
plot(t3t[0:120*60], r3t[0:120*60]),grid()
title("expanded signal")
show()

fL = 1500 # Cutoff frequency 
k = 10 # sinc pulse truncation
tth,ht = sinc_ipol.sinc(fs, fL, k)
plot(tth,ht,"m"),grid()
xlabel("time [sec]"),ylabel("h(t)") 
title("sinc Pulse for Interpolation") 
show()
y3t = convolve(r3t,ht,"same")
plot(t3t[0:120*48],y3t[0:120*48]),grid() 
xlabel("time [sec]"),ylabel("y120t") 
title("Signal 120x Expanded and Interpolated") 
show()

