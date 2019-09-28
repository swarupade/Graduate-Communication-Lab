# Module: comsig.py 
# Class and method definitions for COMmunication Signals 
import numpy as np 
import copy
class sigWave: 
    """ Class for ’waveform’ (CT) signals 
    """ 
    type = "waveform" 
    def __init__(self, sig, Fs=8000, t0=0): 
        """ sig: real or complex-valued waveform samples 
            Fs: sampling rate (default 8000 samples/sec) 
            t0: start time of waveform in seconds (default 0) 
        """ 
        self._sig = np.asanyarray(sig) 
        self._Fs = Fs 
        self._t0 = t0 
        self._shape = np.shape(self._sig) 
        if len(self._shape) > 1: 
            self._Nsamp = len(self._sig[0]) 
        else: 
            self._Nsamp = len(self._sig) 
        self._tlen = self._Nsamp/float(self._Fs) 
        self._tend = self._t0 + (self._Nsamp-1)/float(self._Fs)
   
    # Properties 
    def __len__(self): 
        return self._Nsamp # Returns length in samples 
    def shape(self): 
        return self._shape # Returns shape of signal array 
    def get_Fs(self): 
        return self._Fs # Returns sampling rate 
    def get_t0(self): 
        return self._t0 # Returns start time 
    def get_tlen(self): 
        return self._tlen # Returns length in seconds 
    def get_avgpwr(self): # Returns average power 
        return np.mean(np.power(np.abs(self._sig),2.0)) 
    def get_tend(self): 
        return self._tend # Returns end time 
    def set_t0(self, t0): 
        self._t0 = t0 # Set new start time 
        self._tend = self._t0 + (self._Nsamp-1)/float(self._Fs)
    # Methods 
    def timeAxis(self): # Generate time axis 
        return self._t0 + np.arange(self._Nsamp)/float(self._Fs) 
    def signal(self): # Return the waveform 
        return self._sig 
    def copy(self): # Make a copy of a sigWave object 
        return copy.deepcopy(self)
    def normalized(self): # Normalize the signal to -1,+1 range 
        new_sig = 1.0/np.max(abs(self._sig))*self._sig 
        return sigWave(new_sig, self._Fs, self._t0) 
    def scale(self, factor): # Make a scaled copy of a sigWave object 
        return sigWave(factor*self._sig, self._Fs, self._t0) 
    def pwrx(self, x): # Raise the signal to power x 
        return sigWave(np.power(self._sig, x), self._Fs, self._t0) 
    def apwrx(self, x): # Raise absolute value of signal to power x 
        return sigWave(np.power(np.abs(self._sig), x), self._Fs, self._t0)
    def __add__(self, other): 
        """Add two sigWave signals, sample by sample""" 
        if other == 0: 
            return self 
        assert self._Fs == other._Fs 
        new_t0 = min(self._t0, other._t0) 
        new_Nsamp = 1 + int(round(self._Fs*(max(self._tend,other._tend)-new_t0))) 
        new_tend = new_t0 + (new_Nsamp-1)/float(self._Fs) 
        if self._t0 == new_t0: 
            new_self = np.hstack((self._sig,np.zeros(new_Nsamp-self._Nsamp))) 
            new_other = np.hstack((np.zeros(new_Nsamp-other._Nsamp),other._sig)) 
        else: 
            new_self = np.hstack((np.zeros(new_Nsamp-self._Nsamp),self._sig)) 
            new_other = np.hstack((other._sig,np.zeros(new_Nsamp-other._Nsamp))) 
            new_sig = new_self + new_other 
        return sigWave(new_sig, self._Fs, new_t0) 
    def __or__(self, other): 
        """Concatenate two waveforms""" 
        assert self._Fs == other._Fs 
        new_sig = np.hstack((self._sig,other._sig)) 
        return sigWave(new_sig, self._Fs, self._t0)

class sigSequ: 
    """ Class for ’sequence’ (DT) signals """
    type = "sequence"
    def __init__(self, sig, FB = 100, n0 = 0): 
        """ sig: real or complex-valued sequence values 
        FB: symbol (or Baud) rate (default 100 Baud) 
        n0: start index of sequence (default 0) 
        """ 
        self._sig = np.asanyarray(sig) 
        self._FB = FB 
        self._n0 = n0
    # Properties 
    def __len__(self): 
        return len(self._sig) 
    def get_size(self): 
        return self._sig.size 
    def get_shape(self): 
        return self._sig.shape 
    def get_FB(self): 
        return self._FB 
    def get_n0(self): 
        return self._n0
    # Methods 
    def indexAxis(self): 
        return self._n0 + np.arange(len(self._sig)) 
    def signal(self): 
        return self._sig 
    def copy(self): 
        return copy.deepcopy(self) 
    def scale_offset(self, a, b = 0): 
        """ x[n]_out = a*x[n]_in + b """ 
        return sigSequ(a*self._sig+b, self._FB, self._n0)



