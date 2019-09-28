from matplotlib.pyplot import *
from ascfun import *
from PAM11 import *
import comsig

dn = asc2bin("Test")
an = comsig.sigWave(dn,100,0)
xt = pam11(an,100,10000, "man",[5,0.4])
st = xt.signal()
tt = xt.timeAxis()

matplotlib.pyplot.plot(tt,st)
xlabel("s(t), s(nTb)")