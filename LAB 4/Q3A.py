import numpy as np 
import comsig 
import PAM11 
import showfun1 as showfun
import ascfun 
import matplotlib.pyplot as plt

Fs = 44100 #Sampling rate 
FB = 100 #Baud rate 
N = FB #Number of symbols 
dn = ascfun.asc2bin("Test")
sig_an =comsig.sigSequ(dn,FB,0)
# ***** Set ptype, pparms here *****
ptype= "sinc"
pparms = [10,0.4]
tt , sig_pt = PAM11.pam11(sig_an,FB,Fs,ptype,pparms) #Generate PAM pulse
st = sig_pt.signal()
plt.plot(tt,st)
plt.grid()
plt.title("Polar Binary PAM for string = Test")
plt.show()