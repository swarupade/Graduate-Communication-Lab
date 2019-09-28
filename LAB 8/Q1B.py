# For testing coherent demodulator
import keyfun
import comsig
from pylab import *

an = array([0, 1, 1, 1, 0, 0, 1, 0])     # taking the same data seq as taken in Lab Manual Page(4)
Fb=100
Fs=44100
ptype='rect'
pparms=[]
xtype='coh'
fcparms=[300,-pi/2]

sig_an = comsig.sigSequ(an,Fb)
tt,sig_xt,sig_st = keyfun.askxmtr(sig_an,Fb,Fs,ptype,pparms,xtype,fcparms)

bn, bt,wt,ixn= keyfun.askrcvr(sig_xt,xtype,fcparms,[100,0],ptype,pparms)

print(bn)