import numpy as np
import math
import time
import matplotlib.pyplot as plt
from math import sqrt
from scipy.optimize import fminbound

def cmp_fft (lattice,index):

    fft_n_turns=500

    t = np.arange(fft_n_turns)
    sp = np.fft.fft( [lattice.turns[i].x[index] for i in range (fft_n_turns)] )
    freq = np.fft.fftfreq(t.shape[-1])
    freq_positive=freq[np.where(freq>=0)]
    sp_positive=abs(sp[np.where(freq>=0)])
    peak_frequency_x=freq_positive[sp_positive.argmax()]
    
    original_signal_x=[]
    for i in range(lattice.n_turns):
        original_signal_x.append(lattice.turns[i].x[index])
    original_signal_x=original_signal_x*np.hanning(len(original_signal_x))

    xopt_x,it_x=naff(peak_frequency_x,original_signal_x) 
    #print 'NAFF for horizontal',xopt_x

    sp_y = np.fft.fft( [lattice.turns[i].y[index] for i in range(fft_n_turns)] )
    freq_y = np.fft.fftfreq(t.shape[-1])
    freq_positive_y=freq_y[np.where(freq_y>=0)]
    sp_positive_y=abs(sp_y[np.where(freq_y>=0)])
    peak_frequency_y=freq_positive_y[sp_positive_y.argmax()]
  
    original_signal_y=[]
    for i in range(lattice.n_turns):
      original_signal_y.append(lattice.turns[i].y[index])
    original_signal_y=original_signal_y*np.hanning(len(original_signal_y))

    xopt_y,it_y=naff(peak_frequency_y,original_signal_y)
    #print 'NAFF for vertical',xopt_y
    
    return peak_frequency_x,peak_frequency_y,xopt_x,xopt_y

def cmp_min(f,a,c,b,prec,i):
    i=i+1
    phi = (1 + sqrt(5))/2
    resphi = 2 - phi
    if abs(a - b) < prec:
        return (a + b)/2,i
    d = c + resphi*(b - c)
    if f(d) < f(c):
        return cmp_min(f, c, d, b, prec,i)
    else:
        return cmp_min(f, d, c, a, prec,i)


def expected_value(a,b):
    expected_value=np.sum(np.conj(a[t])*b[t] for t in range (len(a)))/len(a)
    return expected_value

def cmp_exp(freq,N):
    b=[[0]*i for i in range (N)]
    for t in range (N):
        b[t]=np.exp(complex(0,-1)*2*np.pi*freq*t)
    return b

def cmp_sum(freq,N,data):
     a=-abs(expected_value(data,cmp_exp(freq,N)))
     return a


def naff(peak_frequency,data):
     i=0
     y=lambda f:cmp_sum(f,N,data)
     N=len(data)
     small_step=0.005 
     a=peak_frequency-small_step
     b=peak_frequency+small_step
     c=peak_frequency
     return cmp_min(y,a,c,b,1e-8,i)

