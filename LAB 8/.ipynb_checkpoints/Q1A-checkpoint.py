from pylab import *
import comsig
import keyfun
import matplotlib.pyplot as plt

"""
#Coh transmission

an = array([0, 1, 1, 1, 0, 0, 1, 0])                       # taking the same data seq as taken in Lab Manual Page(2)
Fb = 100
Fs=44100
ptype='rect'
pparms=[]
xtype='coh'
fcparms=[300,-pi/2]

sig_an = comsig.sigSequ(an, Fb, 0)
tt,sig_xt,sig_st = keyfun.askxmtr(sig_an,Fb,Fs,ptype,pparms,xtype,fcparms)

plt.plot(tt,sig_xt.signal(), '-r')
plt.xlabel("time [sec]")
plt.ylabel("Amplitude")
plt.title("ASK signal (Coherent): PAM Signal s(t) −−> OOK Signal x(t),\n FB =100 Hz, fc =300 Hz, thetac = -90")
plt.ylim([-2,2])
plt.grid()
plt.show()

"""

# For testing noncoh transmission

an = array([0, 1, 1, 1, 0, 0, 1, 0])     # taking the same data seq as taken in Lab Manual Page(4)
thetacn = 2*pi*rand(len(an))
Fb = 100
Fs=44100
ptype='rect'
pparms=[]
xtype='noncoh'
fcparms=[300,thetacn]

sig_an = comsig.sigSequ(an, Fb, 0)

tt,sig_xt,sig_st = keyfun.askxmtr(sig_an,Fb,Fs,ptype,pparms,xtype,fcparms)
plt.plot(tt,sig_xt.signal(), '-r')
plt.xlabel("time [sec]")
plt.ylabel("Amplitude")
plt.title("ASK signal (Non-Coherent): PAM Signal s(t) −−> OOK Signal x(t),\n FB =100 Hz, fc =300 Hz, Random thetac[n]")
plt.ylim([-2,2])
plt.grid()
plt.show()