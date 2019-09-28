import numpy as np
import wavfun 
import matplotlib.pyplot as plt
import comsig
import showfun1
NTd = 50
dly = 0.5
L= 2
Fb = 32000 
an, Fs = wavfun.wavread("Test_20db_Fb_div_2.wav")
sig_an = comsig.sigWave(an, Fs, 0)
showfun1.showeye(sig_an, Fs, Fb, NTd, [dly, 3, -1.5*L, 1.5*L])