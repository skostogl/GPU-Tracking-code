import sys
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("y.dat")

sp = np.fft.fft( data )
freq = np.fft.fftfreq(len(data))
plt.plot(freq, sp.real, freq, sp.imag)

plt.show()

