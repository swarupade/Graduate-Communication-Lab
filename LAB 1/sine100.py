#sinewave script

from pylab import *
from numpy import *

from ast import literal_eval
Fs = literal_eval(input("Enter sampling rate Fs in Hz: ")) 
fO = 100
tlen = 5e-2
tt = arange(0,round(tlen*Fs))/float(Fs)
st = sin(2*pi*fO*tt)
rt = sign(st)
rdt = diff(hstack((0,rt)))*Fs # Derivative of rt 
 # Integral of rdt
plot(tt[0:248], rdt[0:248])
ylabel("Derivative of rt")
xlabel("time[sec]")
grid()
show()
rdit = cumsum(rdt)/float(Fs)
plot(tt[0:248], rdit[0:248])
ylabel("Integral of derivative of rt")
xlabel("time[sec]")
grid()
show()

            