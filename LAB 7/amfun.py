from pylab import *
import filtfun
import comsig


def amxmtr(tt, sig_mt, xtype, fcparms, fmparms=[], fBparms=[]): 
    """ Amplitude Modulation Transmitter for suppressed (’sc’) 
         and transmitted (’tc’) carrier AM signals 
         >>>>> sig_xt = amxmtr(sig_mt, xtype, fcparms, fmparms, fBparms) <<<<< 
                     where sig_xt: waveform from class sigWave 
                     sig_xt.signal(): transmitted AM signal 
                     sig_xt.timeAxis(): time axis for x(t) 
                     sig_mt: waveform from class sigWave 
                     sig_mt.signal(): modulating (wideband) message signal 
                     sig_mt.timeAxis(): time axis for m(t) 
                     xtype: ’sc’ or ’tc’ (suppressed or transmitted carrier) 
                     fcparms = [fc, thetac] for ’sc’ 
                     fcparms = [fc, thetac, alfa] for ’tc’ 
                     fc: carrier frequency 
                     thetac: carrier phase in deg (0: cos, -90: sin) 
                     alfa: modulation index 0 <= alfa <= 1 
                     fmparms = [fm, km, alfam] LPF at fm parameters 
                               no LPF at fm if fmparms = [] 
                     fm: highest message frequency 
                     km: LPF h(t) truncation to |t| <= km/(2*fm) 
                     alfam: LPF at fm frequency rolloff parameter, linear 
                             rolloff over range 2*alfam*fm
                     fBparms = [fBW, fcB, kB, alfaB] BPF at fcB parameters 
                             no BPF if fBparms = [] 
                     fBW: -6 dB BW of BPF 
                     fcB: center freq of BPF 
                     kB: BPF h(t) truncation to |t| <= kB/fBW 
                     alfaB: BPF frequency rolloff parameter, 
                            linear rolloff over range alfaB*fBW 
        
        """
    mt = sig_mt.signal()
    Fs = int(len(tt)/(tt[-1]-tt[0]))                                        # sampling frequency
    fm = fmparms[0]                                                         # the cutoff frequency for LPF
    k = fmparms[1]                                                          # the truncation parameter of LPF
    alfa = fmparms[2]                                                       # frequency roll-off parameter of LPF
    fc = fcparms[0]                                                         # carrier frequency
    thetac = fcparms[1]                                                      # phase of the carrier signal
    fparms = [fm,0]
    #sig_mt = comsig.sigWave(mt, Fs)
    mt_f, n = filtfun.trapfilt(sig_mt, fparms, k, alfa, 'LPF')                                   # message signal after filtering out with LPF
    if xtype.lower() == 'sc':
        xt = mt_f.signal()*cos(2*pi*fc*tt + thetac)
    else:
        mi = fcparms[2]                                                         # modulation index
        xt = (1+mi*mt_f.signal())*cos(2*pi*fc*tt + thetac)
    return comsig.sigWave(xt,Fs)

def amrcvr(tt, rt_sig, rtype, fcparms, fmparms, fBparms,ftype):
    """
    Amplitude Modulation Receiver for coherent ('coh') reception,
    or absolute value ('abs'), or squaring ('sqr') demodulation,
    or I-Q envelope ('iqabs') detection, or I-Q phase ('iqangle')
    detection.
    >>>>> mthat = amrcvr(tt, rt, rtype, fcparms, fmparms, fBparms) <<<<<
    where mthat: demodulated message signal
    tt: time axis for r(t), mhat(t)
    rt: received AM signal
    rtype: Receiver type from list
    'abs' (absolute value envelope detector),
    'coh' (coherent),
    'iqangle' (I-Q rcvr, angle or phase),
    'iqabs' (I-Q rcvr, absolute value or envelope),
    'sqr' (squaring envelope detector)
    fcparms = [fc, thetac]
    fc: carrier frequency
    thetac: carrier phase in deg (0: cos, -90: sin)
    fmparms = [fm, km, alpham]
    LPF at fm parameters no LPF at fm if fmparms = []
    fm: highest message frequency
    km: LPF h(t) truncation to |t| <= km/(2*fm)
    alpham: LPF at fm frequency rolloff parameter, linear
    rolloff over range 2*alpham*fm
    fBparms = [fBW, fcB, kB, alphaB] BPF at fcB parameters
    no BPF if fBparms = []
    fBW: -6 dB BW of BPF
    fcB: center freq of BPF
    kB: BPF h(t) truncation to |t| <= kB/fBW
    alphaB: BPF frequency rolloff parameter, linear
    rolloff over range alphaB*fBW
    """
    rt = rt_sig.signal()
    Fs = int(len(tt)/(tt[-1]-tt[0]))                                        # sampling frequency
    fm = fmparms[0]                                                         # the cutoff frequency for LPF
    k = fmparms[1]                                                          # the truncation parameter of LPF
    alfa = fmparms[2]                                                       # frequency roll-off parameter of LPF
    fc = fcparms[0]                                                         # carrier frequency
    theta = fcparms[1]                                                      # phase of the carrier signal
    fparms = [fm,0]
    if rtype.lower() == 'coh' :
        vt = 2*rt*cos(2*pi*fc*tt + theta)
        sig_vt = comsig.sigWave(vt, Fs)
        sig_mthat,n = filtfun.trapfilt(sig_vt, fparms, k, alfa,ftype)
    elif rtype.lower() == 'abs':
        vt = abs(rt)
        sig_vt = comsig.sigWave(vt, Fs)
        sig_pt,n = filtfun.trapfilt(sig_vt, fparms, k, alfa, ftype)
        mthat = sig_pt.signal() - mean(sig_pt.signal())
        sig_mthat = comsig.sigWave(mthat, Fs)
    elif rtype.lower() == 'iqangle':
        vit = rt*2*cos(2*pi*fc*tt)
        sig_vit = comsig.sigWave(vit, Fs)
        sig_wit,n = filtfun.trapfilt(sig_vit, fparms, k, alfa, ftype)
        vqt = -rt*2*sin(2*pi*fc*tt)
        sig_vqt = comsig.sigWave(vqt, Fs)
        sig_wqt,n = filtfun.trapfilt(sig_vqt, fparms, k, alfa, ftype)
        mthat = arctan(sig_wqt.signal()/sig_wit.signal())
        sig_mthat = comsig.sigWave(mthat, Fs)
    elif rtype.lower() == 'iqabs':
        vit = rt*2*cos(2*pi*fc*tt)
        sig_vit = comsig.sigWave(vit, Fs)
        sig_wit,n = filtfun.trapfilt(sig_vit, fparms, k, alfa, ftype)
        vqt = -rt*2*sin(2*pi*fc*tt)
        sig_vqt = comsig.sigWave(vqt, Fs)
        sig_wqt,n = filtfun.trapfilt(sig_vqt, fparms, k, alfa, ftype)
        mthat = pow((pow(sig_wit.signal(),2) + pow(sig_wqt.signal(),2)), 0.5)
        sig_mthat = comsig.sigWave(mthat, Fs)
    elif rtype.lower() == 'sqr':
        vt = pow(rt,2)
        sig_vt = comsig.sigWave(vt, Fs)
        sig_wt,n = filtfun.trapfilt(sig_vt, fparms, k, alfa, ftype)
        pt = pow(sig_wt.signal(),2)
        mthat = pt - mean(pt)
        sig_mthat = comsig.sigWave(mthat, Fs)
    else:
        print("rtype is not valid")

    return sig_mthat

def qamxmtr(tt, sig_mt, fcparms, fmparms, ftype):
    """
    Quadrature Amplitude Modulation (QAM) Transmitter with
    complex-valued input/output signals
    >>>>> xt = qamxmtr(tt, mt, fcparms, fmparms) <<<<<
    where xt:complex-valued QAM signal
    tt:time axis for m(t), x(t)
    mt:complex-valued (wideband) message signal
    fcparms = [fc, thetac]
    fc:carrier frequency
    thetac: carrier phase in deg
    fmparms = [fm, km, alpham] for LPF at fm parameters
    fm:highest message frequency (-6dB)
    km:h(t) is truncated to
    |t| <= km/(2*fm) for LPF
    |t| <= km/fBW for BPF
    alpham: frequency rolloff parameter, linear
    rolloff over range
    (1-alpham)*fm <= |f| <= (1+alpham)*fm for LPF
    (1-alpham)*fBW/2 <= |f| <= (1+alpha)*fBW/2 for BPF
    """
    fc, thetaci, thetacq = fcparms[0],fcparms[1],fcparms[2]
    fL, k, alfa = fmparms[0],fmparms[1],fmparms[2]
    fparms = [fc,0]
    sig_xt, n = filtfun.trapfilt_cc(sig_mt, fparms, k, alfa, ftype)
    
    Fs = sig_mt.get_Fs()
    ar = zeros(len(sig_mt.signal()))
    ar = 0j
    xut = zeros(len(sig_mt.signal())) + ar
    xut.real = sig_xt.signal() * cos(2*pi*fc*tt + thetaci)
    xut.imag = sig_xt.signal() * cos(2*pi*fc*tt + thetacq)
    sig_xut = comsig.sigWave(xut, Fs)
    fparms = [fL,0]
    
    return sig_xt

def qamrcvr(tt,sig_rt, fcparms, fmparms,ftype):
    """
    Quadrature Amplitude Modulation (QAM) Receiver with
    complex-valued input/output signals
    >>>>> mthat = qamrcvr(tt, rt, fcparms, fmparms) <<<<<
    where mthat: complex-valued demodulated message signal
    tt:time axis for r(t), mhat(t)
    rt:received QAM signal (real- or complex-valued)
    fcparms = [fc thetac]
    fc:carrier frequency
    thetac: carrier phase in deg
    fmparms = [fm, km, alpham]
    for LPF at fm parameters
    fm:highest message frequency (-6 dB)
    fmparms = [fBW, fBc, km, alpham]
    for BPF at fm parameters
    fBW:BPF -6 dB bandwidth in Hz
    fBc:BPF center frequency (pos/neg) in Hz
    no LPF at fm if fmparms = []
    km:h(t) is truncated to
    |t| <= km/(2*fm) for LPF
    |t| <= km/fBW for BPF
    alpham: frequency rolloff parameter, linear
    rolloff over range
    (1-alpham)*fm <= |f| <= (1+alpham)*fm for LPF
    (1-alpham)*fBW/2 <= |f| <= (1+alpha)*fBW/2 for BPF
    """
    fc, theta = fcparms[0],fcparms[1]
    fL, k, alfa = fmparms[0],fmparms[1],fmparms[2]
    Fs = sig_rt.get_Fs()
    xut = sig_rt.signal()* exp(-1j*(2*pi*fc*tt + theta))
    sig_xut = comsig.sigWave(xut, Fs)
    fparms = [fL,0]
    sig_xt,n = filtfun.trapfilt_cc(sig_xut,fparms, k, alfa, ftype)
    return sig_xt