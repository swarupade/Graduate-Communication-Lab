from pylab import *
import numpy as np
import ascfun
import comsig
import pamfun 
import importlib
importlib.reload(pamfun)

an = ascfun.asc2bin("Test")
an_bipolar = 2*an -1
an_final = hstack((array([0,0]), an_bipolar, array([0,0])))

"""
tt,st=pamfun.pam10(an_final,100,44100,"rect") 
plot(tt,st)
ylim([-2,2])
title('Polar Binary PAM for String \'Test\', FB =100 Baud, Rectangular p(t)' )
grid()
xlabel('time stamp in sec')
ylabel('PAM signal st')
tt,st=pamfun.pam10(an_final,100,44100,"tri") 
plot(tt,st)
ylim([-2,2])
title('Polar Binary PAM for String \'Test\', FB =100 Baud, Triangular p(t)' )
grid()
xlabel('time stamp in sec')
ylabel('PAM signal st')
"""

tt,st=pamfun.pam10(an_final,100,44100,"sinc",[10,4]) 
plot(tt,st)
ylim([-2,2])
title('Polar Binary PAM for String \'Test\', FB =100 Baud, sinc p(t)' )
grid()
xlabel('time stamp in sec')
ylabel('PAM signal st')
