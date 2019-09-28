# File: keyfun.py
# Functions for amplitude/frequency/phase shift keying
# ASK, FSK, PSK and hybrid APSK
from pylab import *
import numpy as np
import comsig
import PAM12
import filtfun
import pamfun
import statistics

def askxmtr(sig_an,Fb,Fs,ptype,pparms,xtype,fcparms):

    """
    Amplitude Shift Keying (ASK) Transmitter for
    Coherent (’coh’) and Non-coherent (’noncoh’) ASK Signals
       >>>>> tt,yyyyyyyyyyyyyyyyyyyysig_xt,sig_st = askxmtr(sig_an,Fs,ptype,pparms,xtype,fcparms)          <<<<<
             where sig_xt: waveform from class sigWave
             sig_xt.signal(): transmitted ASK signal, sampling rate Fs
             x(t) = s(t)*cos(2*pi*fc*t+(pi/180)*thetac)
             sig_xt.timeAxis(): time axis for x(t), starts at t=-TB/2
             sig_st: waveform from class sigWave
             sig_st.signal(): baseband PAM signal s(t) for ’coh’
             sig_st.signal(): st = sit + 1j*sqt for ’noncoh’
             sit: PAM signal of an*cos(pi/180*thetacn)
             sqt: PAM signal of an*sin(pi/180*thetacn)
             xtype: Transmitter type from list {’coh’,’noncoh’}
             sig_an: sequence from class sigSequ
             sig_an.signal() = [an] for {’coh’}
             sig_an.signal() = [[an],[thetacn]] for {’noncoh’}
             an: N-symbol DT input sequence a_n, 0<=n<N
             thetacn: N-symbol DT sequence theta_c[n] in degrees,
             used instead of thetac for {’noncoh’} ASK
             sig_an.get_FB(): baud rate of a_n (and theta_c[n]), TB=1/FB
             Fs: sampling rate of x(t), s(t)
             ptype: pulse type from list
             [’man’,’rcf’,’rect’,’rrcf’,’sinc’,’tri’]
             pparms = [] for {’man’,’rect’,’tri’}
             pparms = [k, alpha] for {’rcf’,’rrcf’}
             pparms = [k, beta] for {’sinc’}
             k: "tail" truncation parameter for {’rcf’,’rrcf’,’sinc’}
             (truncates at -k*TB and k*TB)
             alpha: Rolloff parameter for {’rcf’,’rrcf’}, 0<=alpha<=1
             beta: Kaiser window parameter for {’sinc’}
             fcparms = [fc, thetac] for {’coh’}
             fcparms = [fc] for {’noncoh’}
             fc: carrier frequency in Hz
             thetac: carrier phase in deg (0: cos, -90: sin)
    """
    
    xtype = xtype.lower()
    ptype = ptype.lower()
    if xtype == 'coh':
        fc=fcparms[0]
        tt, sig_st = PAM12.pam12(sig_an, Fb, Fs, ptype, pparms=[])
        sig_st.signal()[np.where(sig_st.signal()<0)] = 0
        thetac = fcparms[1]
        xt = sig_st.signal()*np.cos(2*np.pi*fc*tt + thetac)
        sig_xt = comsig.sigWave(xt, Fs)
    elif xtype == 'noncoh':
        fc = fcparms[0]
        thetacn = fcparms[1]
        it = sig_an.signal()*np.cos(thetacn)
        sig_it = comsig.sigSequ(it, Fb)
        qt = sig_an.signal()*np.sin(thetacn)
        sig_qt = comsig.sigSequ(qt, Fb)
        tt, sig_sit = PAM12.pam12(sig_it, Fb, Fs, ptype, pparms=[])
        tt, sig_sqt = PAM12.pam12(sig_qt, Fb, Fs, ptype, pparms=[])
        st = sig_sit.signal() + 1j*sig_sqt.signal()
        xt = np.real(st*np.exp(1j*2*np.pi*fc*tt))
        sig_xt = comsig.sigWave(xt, Fs)
        sig_st = comsig.sigWave(st, Fs)
    else:
        print("xtype is incorrect")
    
    return tt, sig_xt, sig_st


def askrcvr(tt,sig_rt,Fs,rtype,fcparms,FBparms,ptype,pparms):

    """
        Amplitude Shift Keying (ASK) Receiver for
        Coherent (’coh’) and Non-coherent (’noncoh’) ASK Signals
        >>>>> sig_bn,sig_bt,sig_wt,ixn =
        askrcvr(sig_rt,rtype,fcparms,FBparms,ptype,pparms) <<<<<
        where sig_bn: sequence from class sigSequ
        sig_bn.signal(): received DT sequence b[n]
        sig_bt: waveform from class sigWave
        sig_bt.signal(): received ’CT’ PAM signal b(t)
        sig_wt: waveform from class sigWave
        sig_wt.signal(): wt = wit + 1j*wqt
        wit: in-phase component of b(t)
        wqt: quadrature component of b(t)
        ixn: sampling time indexes for b(t)->b[n], w(t)->w[n]
        sig_rt: waveform from class sigWave
        sig_rt.signal(): received (noisy) ASK signal r(t)
        sig_rt.timeAxis(): time axis for r(t)
        rtype: receiver type from list [’coh’,’noncoh’]
        fcparms = [fc, thetac] for {’coh’}
        fcparms = [fc] for {’noncoh’}
        fc: carrier frequency in Hz
        thetac: carrier phase in deg (0: cos, -90: sin)
        FBparms = [FB, dly]
        FB: baud rate of PAM signal, TB=1/FB
        dly: sampling delay for b(t)->b[n], fraction of TB
        sampling times are t=n*TB+t0 where t0=dly*TB
        ptype: pulse type from list
                [’man’,’rcf’,’rect’,’rrcf’,’sinc’,’tri’]
        pparms = [] for ’man’,’rect’,’tri’
        pparms = [k, alpha] for {’rcf’,’rrcf’}
        pparms = [k, beta] for {’sinc’}
        k: "tail" truncation parameter for {’rcf’,’rrcf’,’sinc’}
           (truncates at -k*TB and k*TB)
        alpha: Rolloff parameter for {’rcf’,’rrcf’}, 0<=alpha<=1
        beta: Kaiser window parameter for {’sinc’}
   """
    
    ptype = ptype.lower()
    rtype = rtype.lower()
    FBparms = FBparms
    pparms = pparms
    if rtype == 'coh':
        fc=fcparms[0]
        thetac = fcparms[1]
        rt = 2*sig_rt.signal()*np.cos(2*np.pi*fc*tt + thetac)
        bn, bt, ixn = pamfun.pamrcvr10(tt, rt, Fs, FBparms, ptype, pparms)
        wt = bt
        bn[np.where(abs(bn)>0.04)] = 1
        bn[np.where(abs(bn)<=0.04)] = 0
      
    
    elif rtype == 'noncoh':
        an = sig_rt.signal()
        fc = fcparms[0]
        thetacn = fcparms[1]
        vt = 2*sig_rt.signal()*np.exp(-1j*2*np.pi*fc*tt)
      
        bni, wit, ixni, = pamfun.pamrcvr10(tt, np.real(vt), Fs, FBparms, ptype, pparms)
        bnq, wqt, ixnq, = pamfun.pamrcvr10(tt, np.imag(vt), Fs, FBparms, ptype, pparms)
        bt = (wit**2 + wqt**2)**0.5
        N = np.ceil(FBparms[0]*(tt[-1]-tt[0]))
        ixn = np.array(np.around((np.arange(N)+0.5+FBparms[1])*Fs/FBparms[0]),int)
        bn = bt[ixn]
        print(2)
        bn[np.where(abs(bn)>0.06)] = 1
        bn[np.where(abs(bn)<=0.06)] = 0
        wt = wit + 1j*wqt
    
    else:
        print("xtype is incorrect")
    bn = np.array(bn,np.int8)
    return bn, bt, wt, ixn

def fskxmtr(M,dnthcn,FB,Fs,ptype,pparms,xtype,fcparms):
    """
    M-ary Frequency Shift Keying (FSK) Transmitter for
    Choherent ('coh'), Non-coherent ('noncoh'), and
    Continuous Phase ('cpfsk') FSK Signals
    >>>>> tt,xt = fskxmtr(M,dnthcn,FB,Fs,ptype,pparms,xtype,fcparms) <<<<<
    where  tt:      time axis for x(t), starts at t=-TB/2
    xt:      transmitted FSK signal, sampling rate Fs
    M:       number of distinct symbol values in d[n]
    xtype:   Transmitter type from set {'coh','noncoh','cpfsk'}
    dnthcn = [dn]           for ['coh','cpfsk']
    dnthcn = [dn, thetacn]  for ['noncoh']
    dn:      M-ary (0,1,..,M-1) N-symbol DT input sequence d_n
    thetacn: N-symbol DT sequence theta_c[n] in degrees,
    used instead of thetac0..thetacM-1 for {'noncoh'} FSK
    FB:      baud rate of d_n (and theta_c[n]), TB=1/FB
    Fs:      sampling rate of x(t)
    ptype:   pulse type from set
    {'man','rcf','rect','rrcf','sinc','tri'}
    pparms = []         for {'man','rect','tri'}
    pparms = [k alpha]  for {'rcf','rrcf'}
    pparms = [k beta]   for {'sinc'}
    k:       "tail" truncation parameter for {'rcf','rrcf','sinc'}
    (truncates p(t) to -k*TB <= t < k*TB)
    alpha:   Rolloff parameter for {'rcf','rrcf'}, 0<=alpha<=1
    beta:    Kaiser window parameter for {'sinc'}
    fcparms = [[fc0,fc1,...,fcM-1],[thetac0,thetac1,...,thetacM-1]]
    for {'coh'}
    fcparms = [fc0,fc1,...,fcM-1]  for {'noncoh'}
    fcparms = [fc, deltaf]         for {'cpfsk'}
    fc0,fc1,...,fcM-1:   FSK (carrier) frequencies
    for {'coh','noncoh'}
    thetac0,thetac1,...,thetacM-1: FSK (carrier) phases in deg
    (0: cos, -90: sin) for {'coh'}
    fc:      carrier frequency for {'cpfsk'}
    deltaf:  frequency spacing for {'cpfsk'}
    for dn=0 -> fc, dn=1 -> fc+deltaf,
    dn=2 -> fc+2*deltaf, etc
    """
    if xtype == 'cpfsk':
        fc0 = fcparms[0]
        fc1 = fcparms[1]
    elif xtype == 'coh':
        fc0 = fcparms[0]
        fc1 = fcparms[1]
    else:
        fc0 = fcparms[0][0]
        fc1 = fcparms[0][1]

    if xtype == 'coh':
        di = dnthcn[0]      
    else:
        di = dnthcn          
    N = len(di)                
    TB = 1/float(FB)           
    ixL = ceil(-Fs*0.5*TB)     
    ixR = ceil(Fs*(N-0.5)*TB)  
    tt = arange(ixL,ixR)/float(Fs)  
    
    # ***** Conversion from DT a_n to CT a_s(t) *****
    ast = zeros(len(tt))      
    ast1 = zeros(len(tt))      
    ix = array(around(Fs*arange(0,N)*TB),int)
    
    
    ast[ix-int(ixL)] = Fs*di   
    di1 = array([0 if k>0 else 1 for k in di]) 
    ast1[ix-int(ixL)] = Fs*di1
    if xtype == 'noncoh':
        thetac_n = [float((fcparms[1][k])*pi/180) for k in range(0,len(fcparms[1]))]  
        print(thetac_n)
    
    # ***** Set up PAM pulse p(t) *****
    ptype = ptype.lower()     
    
    # Set left/right limits for p(t)
    if (ptype=='rect'):
        kL = -0.5; kR = -kL
    else:
        kL = -1.0; kR = -kL    
    ixpL = ceil(Fs*kL*TB)      
    ixpR = ceil(Fs*kR*TB)      
    ttp = arange(ixpL,ixpR)/float(Fs)  
    pt = zeros(len(ttp))       
    if (ptype=='rect'):       
        ix = where(logical_and(ttp>=kL*TB, ttp<kR*TB))[0]
        pt[ix] = ones(len(ix))
    elif(ptype == 'tri'):
        pt = [(1+((i*1)/TB)) if (1+((i*1)/TB))<1 else (1-((i*1)/TB)) for i in list(ttp)]
    elif(ptype == 'sinc'):
        k=pparms[0]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = zeros(len(ttp))
        beta=pparms[1]
        pt = [(sin(pi*q/TB)/(pi*q/TB)) if q!=0 else 1 for q in list(ttp)]
        pt = pt*kaiser(len(pt),beta)
    elif(ptype=='rcf'):
        k = pparms[0]
        alpha = pparms[1]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = zeros(len(ttp))
        j = 0
        pt = [(sin((pi*q)/TB)/((pi*q)/TB)*(cos(pi*alpha*q/TB)/(1-pow(2*alpha*q/TB,2)))) if q!=0 else 1 for q in list(ttp)]
        ix_nz = where(ttp==-(TB/(2*alpha)))[0]
        ix_pz = where(ttp==(TB/(2*alpha)))[0]
        pt[ix_nz] = (pt[ix_nz+1] + pt[ix_nz-1])/2
        pt[ix_pz] = (pt[ix_pz+1] + pt[ix_pz-1])/2
        
    elif(ptype=='rrcf'):
        k = pparms[0]
        alpha = pparms[1]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = zeros(len(ttp))
        ixp = where(ttp!=0)[0]
        p = pi*ttp[ixp]/TB
        a = alpha*ttp[ixp]/TB
        pt[ixp] = ((sin((1-alpha)*p) + 4*a*cos((1+alpha)*p))/((1-np.power(4*a,2))*pi*ttp[ixp]))*(TB)      
        ix1 = where(ttp==0)[0]
        pt[ix1] = (1-alpha+(4*alpha/pi))
        ix2 = where(ttp==TB/(4*alpha))
        ix3 = where(ttp==-TB/(4*alpha))
        if ix2:    
            pt[ix2] = (alpha/2**0.5)*(((1+2/pi)*sin(pi/(4*alpha))) + ((1-2/pi)*cos(pi/(4*alpha))))
            pt[ix3] = (alpha/2**0.5)*(((1+2/pi)*sin(pi/(4*alpha))) + ((1-2/pi)*cos(pi/(4*alpha))))
    elif(ptype=='man'):
        kL = -0.5; kR = -kL
        ixpL = ceil(Fs*kL*TB)     
        ixpR = ceil(Fs*kR*TB)     
        ttp = arange(ixpL,ixpR)/float(Fs) 
        pt = zeros(len(ttp))
        ix = where(logical_and(ttp>=kL*TB, ttp<0))[0]
        pt[ix] = -(ones(len(ix)))
        ix = where(logical_and(ttp>=0, ttp<kR*TB))[0]
        pt[ix] = (ones(len(ix)))
    else:
        print("ptype '%s' is not recognized" % ptype)
    st = convolve(ast,pt)/float(Fs)  
    st = st[-int(ixpL):int(ixR-ixL-ixpL)]
    
    st1 = convolve(ast1,pt)/float(Fs)  
    st1 = st1[-int(ixpL):int(ixR-ixL-ixpL)]
    yt = array([])
    if xtype == 'noncoh':
        tt1 = tt[0:int(TB*Fs)]
        for i in range(0,len(di)):
            dt_sym = st[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc0*tt1)+thetac_n[i]) + st1[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc1*tt1)+thetac_n[i]) #data for each symol
            yt = append(yt,dt_sym)
    elif xtype == 'coh':
        tt1 = tt[0:int(TB*Fs)]
        for i in range(0,len(di)):
            dt_sym = st[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc0*tt1)+float(dnthcn[1]*pi/180)) + st1[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc1*tt1)+float(dnthcn[1]*pi/180)) #data for each symol
            yt = append(yt,dt_sym)
    elif xtype == 'cpfsk':
        tt1 = tt[int(TB*Fs/2):int(3*TB*Fs/2)]
        for i in range(0,len(di)):
            dt_sym = st[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc0*tt1)+float(fcparms[2]*pi/180)) + st1[i*(len(tt1)):(i+1)*len(tt1)]*cos((2*pi*fc1*tt1)+float(fcparms[2]*pi/180)) #data for each symol
            yt = append(yt,dt_sym)
    return tt,yt


def fskrcvr(tt, rt, FBparms, FCparms,rtype, ptype, pparms=[]):
    """
    Pulse amplitude modulation receiver with matched filter:
    r(t) -> b(t) -> bn.
    V1.0 for 'man', 'rcf', 'rect', 'rrcf', 'sinc', and 'tri'
    pulse types.
    >>>>> bn, bt, ixn = pamrcvr10(tt, rt, FBparms, ptype, pparms) <<<<<
    where  tt:    time axis for r(t)
    rt:    received (noisy) PAM signal r(t)
    FBparms: = [FB, dly]
    FB:    Baud rate of PAM signal, TB=1/FB
    dly:   sampling delay for b(t) -> b_n as a fraction of TB
    sampling times are t=n*TB+t0 where t0 = dly*TB
    ptype: pulse type from list
    ('man','rcf','rect','rrcf','sinc','tri')
    pparms not used for 'man','rect','tri'
    pparms = [k, alpha]  for 'rcf','rrcf'
    pparms = [k, beta]  for 'sinc'
    k:     "tail" truncation parameter for 'rcf','rrcf','sinc'
    (truncates p(t) to -k*TB <= t < k*TB)
    alpha: rolloff parameter for ('rcf','rrcf'), 0<=alpha<=1
    beta:  Kaiser window parameter for 'sinc'
    bn:    received DT sequence after sampling at t=n*TB+t0
    bt:    received PAM signal b(t) at output of matched filter
    ixn:   indexes where b(t) is sampled to obtain b_n
    """
    Fs = (len(tt)-1)/(tt[-1]-tt[0])
    thetac = 90
    theta = float(thetac*pi/180)
    if rtype == 'coh':
        Fc = FCparms[0]
    elif rtype == 'noncoh' :
        Fc = FCparms[0][0]
    elif rtype == 'cpfsk':
        Fc = FCparms[0]
    else:
        print('Unrecognized rtype')
    carr_sig = cos((2*pi*Fc*tt)+theta)
    rt = rt * carr_sig
    FB = FBparms[0]
    t0 = FBparms[1]
    
    # ***** Set up matched filter response h_R(t) *****
    TB = 1/float(FB)           
    ptype = ptype.lower()      
    
    # Set left/right limits for p(t)
    if (ptype=='rect'):        
        kL = -0.5; kR = -kL
    else:
        kL = -1.0; kR = -kL    
    ixpL = ceil(Fs*kL*TB)      
    ixpR = ceil(Fs*kR*TB)     
    ttp = arange(ixpL,ixpR)/float(Fs)  
    pt = zeros(len(ttp))       
    if (ptype=='rect'):        
        ix = where(logical_and(ttp>=kL*TB, ttp<kR*TB))[0]
        pt[ix] = ones(len(ix))
    elif (ptype=='tri'):
        pt = [(1+((i*1)/TB)) if (1+((i*1)/TB))<1 else (1-((i*1)/TB)) for i in list(ttp)]
    elif (ptype=='rcf'):
        k = pparms[0]
        alpha = pparms[1]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = [(sin((pi*q)/TB)/((pi*q)/TB)*(cos(pi*alpha*q/TB)/(1-pow(2*alpha*q/TB,2)))) if q!=0 else 1 for q in list(ttp)]
        ix_nz = where(ttp==-(TB/(2*alpha)))[0]
        ix_pz = where(ttp==(TB/(2*alpha)))[0]
        if (ix_nz !=0):
            pt[ix_nz] = (pt[ix_nz+1] + pt[ix_nz-1])/2
            pt[ix_pz] = (pt[ix_pz+1] + pt[ix_pz-1])/2
    elif (ptype=='sinc'):
        k=pparms[0]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = zeros(len(ttp))
        beta=pparms[1]
        pt = [(sin(pi*q/TB)/(pi*q/TB)) if q!=0 else 1 for q in list(ttp)]
        pt = pt*kaiser(len(pt),beta)
    elif (ptype=='man'):
        kL = -0.5; kR = -kL
        ixpL = ceil(Fs*kL*TB)     
        ixpR = ceil(Fs*kR*TB)      
        ttp = arange(ixpL,ixpR)/float(Fs) 
        pt = zeros(len(ttp))
        ix = where(logical_and(ttp>=kL*TB, ttp<0))[0]
        pt[ix] = -(ones(len(ix)))
        ix = where(logical_and(ttp>=0, ttp<kR*TB))[0]
        pt[ix] = (ones(len(ix)))
    elif (ptype=='rrcf'):
        k = pparms[0]
        alpha = pparms[1]
        kL = -1.0*k
        kR = 1.0*k
        ixpR = ceil(Fs*kR*TB)
        ixpL = ceil(Fs*kL*TB)
        ttp = arange(ixpL,ixpR)/float(Fs)
        pt = zeros(len(ttp))
        ixp = where(ttp!=0)[0]
        p = pi*ttp[ixp]/TB
        a = alpha*ttp[ixp]/TB
        pt[ixp] = ((sin((1-alpha)*p) + 4*a*cos((1+alpha)*p))/((1-np.power(4*a,2))*pi*ttp[ixp]))*(TB)      
        ix1 = where(ttp==0)[0]          
        pt[ix1] = (1-alpha+(4*alpha/pi))
        ix2 = where(ttp==TB/(4*alpha))[0]
    else:
        print("ptype '%s' is not recognized" % ptype)
    hRt = pt[::-1]            
    hRt = Fs/sum(np.power(pt,2.0))*hRt  

    # ***** Filter r(t) with matched filter h_R(t)
    bt = convolve(rt,hRt)/float(Fs)  
    bt = bt[int(-ixpL):int(len(tt)-ixpL)]  
    
    # ***** Sample b(t) at t=n*TB+t0 to obtain b_n *****
    N = ceil(FB*(tt[-1]-tt[0])) 
    ixn = array(around((arange(N)+0.5+t0)*Fs/FB),int)
    
    ix = where(logical_and(ixn>=0,ixn<len(tt)))[0]
    ixn = ixn[ix]              
    bn = bt[ixn] 
    bn = [1 if abs(k)>0.03 else 0 for k in bn]
    return bn, bt, ixn
   