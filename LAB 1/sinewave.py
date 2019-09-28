#sinewave module
from pylab import *
import pcmfun as pcm
from ast import literal_eval

Fs = literal_eval(input("Enter sampling rate Fs in Hz: "))
fm = 100
tlen = 5e-2

tt = arange(0, round(tlen*Fs))/float(Fs)
xt = sin(2*pi*fm*tt)
title("original message signal")
plot(tt, xt), grid()
show()

quantized_mt, mt_code = pcm.mt2pcm(xt, 8)
plot(tt[0:200],mt_code[0:200]), grid()
title("Quantized bit string")
show()

mt = pcm.pcm2mt(mt_code, 8)
plot(tt[0:400],mt[0:400]), grid()
title("Converted Analog signal from the PCM")
show()
