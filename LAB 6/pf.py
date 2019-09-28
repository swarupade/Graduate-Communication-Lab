# this module will be imported in the into your flowgraph

# File: ptfun
# Functions for gnuradio-companion PAM p(t) generation
import numpy as np
from pylab import *
def pampt(sps, ptype, pparms=[], plot='', duty=1):
    """
    PAM pulse p(t) = p(n*TB/sps) generation
    >>>>> pt = pampt(sps, ptype, pparms) <<<<<
    where  sps: samples per symbol (Fs/FB)
        ptype: pulse type ('rect', 'sinc', 'tri', 'man', 'rcf', 'rrcf')
        pparms not used for 'rect', 'tri'
        pparms = [k, beta] for sinc
            k:     "tail" truncation parameter for (truncates p(t) to -k*sps <= n < k*sps)
            beta:  Kaiser window parameter for 'sinc'
        pt: pulse p(t) at t=n*TB/sps
    Note: In terms of sampling rate Fs and baud rate FB, sps = Fs/FB
    """
    if ptype is 'rect':
        pt = np.ones(sps)
    elif ptype is 'tri':
        triarray = np.arange(0,1,(1/float(sps)))[1:]
        pt = np.concatenate([triarray,[1],triarray[::-1]])
    elif ptype is 'sinc':
        k = pparms[0]
        beta = pparms[1]
        nn = np.arange(-k*sps,k*sps) # was (-2*k*sps,2*k*sps)
        pt = sinc((1/float(sps))*nn)
        pt = pt*kaiser(len(pt),beta)
    elif ptype is 'man':
        if(sps % 2 == 0): # is even....
            pt = concatenate([-1*ones(int(sps/2)),ones(int(sps/2))])
        else: # is odd....
            pt = concatenate([-1*ones(int(floor(sps/2))),[0],ones(int(floor(sps/2)))])
    elif ptype is 'rcf':
        k = pparms[0]
        alpha = pparms[1]
        nn = np.arange(-k*sps,k*sps)
        tt = nn/float(sps)
        pt=[]
        for t in tt:
            rcft_num = sin(pi*t)*cos(pi*alpha*t)
            rcft_den = (pi*t)*(1-pow(2*alpha*t,2))
            if (rcft_den == 0.0):
                rcft_num=pt[-1]
                rcft_den=1
            rcft = divide(rcft_num,float(rcft_den))
            pt = concatenate([pt,[rcft]])
    elif ptype is 'rrcf':
        k = pparms[0]
        alpha = pparms[1]
        nn = np.arange(-k*sps,k*sps)
        tt = nn/float(sps)
        pt=[]
        for t in tt:
            if(t==0):
                rcft_num = (1-alpha+(4*alpha/pi))
                rcft_den = 1
            elif(abs(t)==1/float(4*alpha)):
                rcft_num = alpha*((1+2/pi)*sin(pi/float(4*alpha))+(1-2/pi)*cos(pi/float(4*alpha)))
                rcft_den = pow(2,0.5)
            else:
                rcft_num = (sin(pi*t*(1-alpha))+(4*alpha*t)*cos(pi*t*(1+alpha)))
                rcft_den = pi*t*(1-pow(4*alpha*t,2))
            rcft = divide(rcft_num,float(rcft_den))
            pt = concatenate([pt,[rcft]])
    else:
        print("ERROR: ptype '",ptype,"' not recognized")
        return 0

    if(duty!=1):
        if(ptype=='tri'):
            sps = sps*2
        elif(ptype=='rcf' or ptype=='sinc' or ptype=='rrcf'):
            sps = sps/duty
        widthbuff = zeros(int(((sps/float(duty))-len(pt))/float(2)))
        pt = concatenate([widthbuff,pt,widthbuff])

    return(pt)

def pamhRt(sps, ptype, pparms=[]):
    """
        PAM normalized matched filter (MF) receiver filter
        h_R(t) = h_R(n*TB/sps) generation
        >>>>> hRt = pamhRt(sps, ptype, pparms) <<<<<
        where sps:
            ptype: pulse type from list
                ('man', 'rcf', 'rect', 'rrcf', 'sinc', 'tri')
            pparms:
                pparms not used for 'man', 'rect', 'tri'
                pparms = [k, alpha] for 'rcf', 'rrcf'
                pparms = [k, beta] for 'sinc'
                k: "tail" truncation parameter for 'rcf', 'rrcf', 'sinc'
                    (truncates p(t) to -k*sps <= n < k*sps)
                alpha: Rolloff parameter for 'rcf', 'rrcf', 0 <= alpha <= 1
                beta: Kaiser window parameter for 'sinc'
            hRt: MF impulse response h_R(t) at t=n*TB/sps
        Note: In terms of sampling rate Fs and baud rate FB, sps = Fs/FB
    """
    pt = pampt(int(sps), ptype, pparms)
    hrt = multiply(pt,1/float(np.sum(np.power(pt,2))))
    hrt = hrt[::-1]
    return hrt
