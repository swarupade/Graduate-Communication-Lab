import numpy as np
import ascfun
import wavfun

rt,fs= wavfun.wavread('pr1sig401.wav')
fb = 100
tb = 1/fb
nsb = float(fs)*tb         # no. of samples in one bit
nbr = round(len(rt)/nsb)      # no. of bits received
nbr = nbr - nbr%8          # modifying number of bits received to be multiple of 8
bn = np.array([])
st = 2*rt     #  non-normalizing the signal
for i in np.arange(nbr):
    bn = np.append(bn,round(st[int(nsb*(0.5 + i))]))   
dn = (bn/2 + 1)%2         # Generating the encoded sequence in zeros and ones using the concept in lab manual
txt = ascfun.bin2asc(dn,bits = 8)   # Decoding the text from the dn sequence
print(txt)