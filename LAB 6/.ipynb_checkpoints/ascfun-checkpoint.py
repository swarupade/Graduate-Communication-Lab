from pylab import *
import numpy as np
def asc2bin(txt, bits=8): 
    """ 
    ASCII text to serial binary conversion 
    >>>>> dn = asc2bin(txt, bits) <<<<< 
    where txt         input text string 
          abs(bits)         bits per char, default=8 
          bits > 0          LSB first parallel to serial 
          bits < 0          MSB first parallel to serial dn binary output sequence 
    """
    txtnum = array([ord(c) for c in txt]) # int array 
    if bits > 0: # Neg powers of 2, increasing exp  
        p2 = np.power(2.0,arange(0,-bits,-1)) 
    else: # Neg powers of 2, decreasing exp 
        p2 = np.power(2.0,1+arange(bits,0)) 
    B = array(mod(array(floor(outer(txtnum,p2)),int),2),int8) 
         # Rows of B are bits of chars 
    dn = reshape(B,B.size) 
    return dn       # Serial binary output

def bin2asc(dn, bits, flg=1):
    
    n = int(len(dn)) # number of bits
    N = int(floor(n/float(abs(bits)))) # number of letters

    bitString = []
    for i in range(N):
        bitString.append(dn[(i*abs(bits)):((i+1)*abs(bits))])
    if bits < 0:
        for i in range(N):
            bitString[i] = bitString[i][::-1]
    ASCIIString = []
    for i in range(N):
        value=0
        for j in range(abs(bits)):
            value = value + (int(bitString[i][j])<<j)
        ASCIIString.append(value)
    textString=""
    for i in range(N):
        textString=textString+chr(ASCIIString[i])

    return(textString,ASCIIString)
