from pylab import *
def sinc(Fs, fL, k):
    # Time axis 
    ixk = int(round(Fs*k/float(2*fL)))
    tth = arange(-ixk,ixk+1)/float(Fs)
    # Sinc pulse
    ht = 2.0*fL*ones(len(tth)) 
    ixh = where(tth != 0.0)[0] 
    ht[ixh] = sin(2*pi*fL*tth[ixh])/(pi*tth[ixh]) 
    return tth,ht
    if __name__ == "__main__": 
        sinc()
