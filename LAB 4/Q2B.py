import numpy as np
import wavfun 
import matplotlib.pyplot as plt
import comsig
import showfun1
NTd = 50
dly = 0.5
L=4
#FB = [300,350,400] for pamsig401
#FB = [100, 150, 160, 180, 200, 220] for pamsig402
Fs = 32000   # for Q2C pamsig403
FB = 100     # for Q2C pamsig403
an, framerate = wavfun.wavread("pamsig403_rcvrhighFl.wav")
tt=np.arange(len(an))/float(framerate)
#plt.plot(tt[:5000],an[0:5000])
#plt.xlabel("Amplitude")
#plt.ylabel("time")
#plt.grid()
#plt.show()

xt = comsig.sigWave(an, Fs)
showfun1.showeye(xt, Fs, FB, NTd, [dly, 3, -0.25*L, 0.25*L])
#showfun1.showft(tt,xt,Fs,[-2000,2000,-60])

