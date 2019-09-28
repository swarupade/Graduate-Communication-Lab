import wavfun
import numpy as np
import showfun1
import comsig

xt, Fs = wavfun.wavread("pr1sig401.wav")
N = Fs
sig_xt = comsig.sigWave(xt, Fs, 0)
showfun1.showpsd(sig_xt, Fs, [-200,200,-60],N)