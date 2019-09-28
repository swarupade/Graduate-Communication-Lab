# File: ftpam01.py 
# Script file that accepts an ASCII text string as input and 
# produces a corresponding binary unipolar flat-top PAM signal 
# s(t) with bit rate Fb and sampling rate Fs. 
# The ASCII text string uses 8-bit encoding and LSB-first 
# conversion to a bitstream dn. At every index value # n=0,1,2,..., dn is either 0 or 1. To convert dn to a 
# flat-top PAM CT signal approximation s(t), the formula 
# s(t) = dn, (n-1/2)*Tb <= t < (n+1/2)*Tb, 
# is used. 

from pylab import * 
import ascfun as af
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal 
import wavfun as wf

Fs = 10000       # Sampling rate for s(t) 
Fb = 100         # Bit rate for dn sequence 
txt = "Test"     # Input text string
bits = 8         # Number of bits 
Tb = 1/float(Fb) # Time per bit 


dn = af.asc2bin(txt, bits)      # >> Convert txt to bitstream dn here << 
print(dn)

N = len(dn)                     # Total number of bits
ixl = round(-0.5*Tb*Fs) 	    # the starting left index of the st
ixr = round((N-0.5)*Fs*Tb)	    # the right most index of the st
tt = arange(ixl,ixr)/float(Fs)  # time axis for st

sdt = diff(hstack((0,dn)))*Fs   # taking differential of dn
sdpt=[]
for i in sdt:                   # adding zeros in between two dn bits 
	sdpt = sdpt + [i] + list(zeros(round(Tb*Fs)-1)) 


sdpt_array = array(sdpt)   	# changing list into numpy array
sdit = cumsum(sdpt_array)/float(Fs)  # integrating the differential of sdt to get sdit

plot(tt[0:14100],sdit[0:14100]), grid()
xlabel("time")
ylabel("s(t)")
title("Flat top PAM signal generated")
show()
wf.wavwrite("Test.wav",Fs,0.999*sdit/float(max(abs(sdit))))
