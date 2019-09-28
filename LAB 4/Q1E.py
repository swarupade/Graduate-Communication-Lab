import numpy as np 
import comsig 
import PAM11 
import showfun1 as showfun
import ascfun
 
Fs = 44100 #Sampling rate 
FB = 100 #Baud rate 
N = FB #Number of symbols 
dn = ascfun.asc2bin("Test")
sig_an =comsig.sigSequ(dn,FB,0)
# ***** Set ptype, pparms here *****
ptype= "sinc"
pparms = [10,0.4]
sig_pt = PAM11.pam11(sig_an,FB,Fs,ptype,pparms) #Generate PAM pulse
st = sig_pt.signal()
tt = sig_pt.timeAxis()
sig_pt.set_t0(sig_pt.get_t0()-(N/2)/float(FB)) #Place center of pulse at t=0 
# ***** Set ff_parms here ***** 
ff_parms = [-200,200,-60]
showfun.showft(tt,sig_pt,Fs,ff_parms) 
#Plot FT of pulse