from pylab import *
import ascfun as af
import wavfun as wf
import numpy as np
import showfun
import importlib
importlib.reload(showfun)

fs, rt  = wf.wavread("pamsig203.wav")
bits = 8
 	
dn=[]
comp_val = (np.max(rt) + np.min(rt))/2.0   # to decide whether the bit is 1 or zero

# making a list of all bits
for j in rt:
	if j > comp_val:
		dn = dn + [1]
	else:
		dn = dn + [0]

# counting the continuous 1s or 0s
k=dn[0] 
n=0  # track the count
lst=[]  # this list will contain the count of continuous 1s or 0s
for i in dn:
	if k == i:
		n=n+1
		continue
	else:
		k=i
		lst = lst + [n]
		n=0

fb=fs/min(lst)


# getting the minimum value of the list and printing it
print("The bit rate is below ")
print(fb)
n=len(rt)      # no. of samples in rt
tlen=(1/fs)*n      # the time interval
tt = arange(0,tlen,1/fs)         #defining time axis
plot(tt,rt)
xlabel('time axis in sec', fontsize=14)
ylabel('amplitude', fontsize=14)
ylim([-2,2])
xlim([0,0.1])
title('The time domain plot of pamsig203.wav file', fontsize=16, color='b')
grid()
showfun.showft(tt,rt,fs,[-500,500,-40])

