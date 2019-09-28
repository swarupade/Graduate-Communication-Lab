import numpy as np 
import comsig 
import PAM11 
import showfun1 as showfun

Fs = 44100 # Sampling rate 
FB = 100 # Baud rate
tlen = 2
NTd = 50 # Number of traces displayed 
N = NTd+10 # Number of data symbols 
L = 2 # Number of data levels 
dly = 0.5 # Trigger delay TB/2 

#####Creating the random polar binary PAM signal of s(t) of length 2 sec and triangular p(t)########

dn = np.around(np.random.rand(tlen*FB)) # Unipolar L-valued random data 
an = (2*dn-float(L/2)) # Polar L-valued DT sequence 
sig_an = comsig.sigSequ(an, FB, 0) # DT sequence class 
tt, sig_st = PAM11.pam11(sig_an, FB, Fs, "tri") # PAM signal, ’tri’ p(t) 

### generating a bandlimited gaussian noise n(t)###             
nn = np.random.randn(np.round(tlen *2 * FB)) # Gaussian noise n(t) 
sig_nn = comsig.sigSequ(nn, 2 *FB, 0)
N = len(sig_nn.signal())
tt, sig_nt = PAM11.pam11(sig_nn, 2 * FB , Fs,'rcf', [20, 0.2])
              
#### Adding them together #####
A = 0.05
rt = sig_st.signal() + (A*sig_nt.signal())
sig_rt = comsig.sigWave(rt, FB, 0)
showfun.showeye(sig_rt, Fs, FB, NTd, [dly, 3, -1.5*L, 1.5*L])
