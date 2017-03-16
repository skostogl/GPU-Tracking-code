import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy.signal import hilbert, chirp

def exp_fit (x,y):
  x = np.array(x, dtype=float) 
  y = np.array(y, dtype=float)

  analytic_signal = hilbert(y)
  amplitude_envelope = np.abs(analytic_signal)
  amplitude_envelope2 = -np.abs(analytic_signal)

  def func(x, a, c, d):
      return a*np.exp(-c*x)+d

  def func2(x, a, c, d):
      return a*np.exp(c*x)+d

  popt, pcov = curve_fit(func, x, amplitude_envelope, p0=(1,1e-6,1))
  popt2, pcov = curve_fit(func2, x, amplitude_envelope2, p0=(1,1e-6,1))

  print "TOP envelope: a = %s , c = %s, d = %s" % (popt[0], popt[1], popt[2])
  print "Bottom envelope: a = %s , c = %s, d = %s" % (popt2[0], popt2[1], popt2[2])
  a = max(y)
  b = min(y)

  f0 =popt[0]+popt[2]
  for i in range (len(x)):
    fi =popt[0]*np.exp(-popt[1]*i)+popt[2]
    y[i] *= f0/fi
  return x,y
