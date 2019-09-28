import numpy as np 
import comsig 
import showfun 
import importlib
importlib.reload(showfun)

Fs = 44100 # Sampling rate 
fa, fb = 140, 164 # Frequencies fa, fb 
tlen = 0.25 # Length of t-axis in sec 
tt = np.arange(0,round(Fs*tlen))/float(Fs) # Time axis 
xt = np.sin(2*np.pi*fa*tt)+0.01*np.cos(2*np.pi*fb*tt) # Linear combination of sinusoids 
sig_xt = comsig.sigWave(xt, Fs, 0) # Waveform from class sigWave 
showfun.showft(sig_xt,[-200,200,-60]) # Display X(f), using ff_lim