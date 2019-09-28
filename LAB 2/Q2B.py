import numpy as np
import showfun
import importlib
importlib.reload(showfun)
import random
import pamfun
import matplotlib.pyplot as plt
Fs = 44100
ptype = 'rect'
pparms =[10,4]
FB = 100
dn = np.random.rand(int(np.round(FB/2))) # Random sequence, uniform in [0...1) 
dn_unipolar = np.array(np.floor(2*dn),int) # Random unipolar binary sequence in {0,1} 
sig_an = 2*dn_unipolar-1                           # Random polar binary sequence in {-1,+1}
tt,st= pamfun.pam10(sig_an,FB,Fs,ptype,pparms)   
plt.plot(tt,st)
plt.title('Polar Binary PAM for random discrete time sequence, FB =%s Baud, rect p(t)' %FB)
plt.grid()
plt.xlabel('time stamp in sec')
plt.ylabel('PAM signal st')
showfun.showft(tt,st,Fs,[-1000,1000,-60])