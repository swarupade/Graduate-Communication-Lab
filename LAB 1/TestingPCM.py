#pam and pcm combined for Q3 part b
from pylab import *
import ascfun as af
import wavfun as wf
import pcmfun as pfun

fs, rt = wf.wavread("pcm_test01.wav")

fb = 24000
tb = 1/float(fb)

fO = 233.3
tlen = 5e-2
tt = arange(0,round(tlen*fs))/float(fs)
st = sin(2*pi*fO*tt)
plot(tt, st), grid()
title("Ideal Waveform")
show()

bits = 8
n = int(floor(len(rt)/float(fs)/tb)) 	# number of received bits

rt = list(rt)  	# changing rt into list type

######### getting sample of rt signal #######

dnhat=[]
for i in range(n):
	d_prime = rt[i*round(fs*tb):(i+1)*round(fs*tb)]
	avg = sum(d_prime) / round(fs*tb)     # averaging out the one bit window and the comapring
	if avg > 0.5:
		dnhat = dnhat + [1]
	else:
		dnhat = dnhat + [0]

#####################################################
dnhat = array(dnhat,int8)	 # converting list into binary array
rt = pfun.pcm2mt(dnhat,8)

plot(tt[0:200], rt[0:200]), grid()
title("Converted Analog signal from PCM")
show()

#wf.wavwrite("pcm_sig02_analog.wav",fs,0.999*rt/float(max(abs(rt))))

