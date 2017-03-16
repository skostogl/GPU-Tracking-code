import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from modules.NAFF import *
from random import seed
from random import gauss
from random import random

t = np.arange(0.0,1.0,0.001)
data=[]

########## Noise #######
#seed(1)
#series = [5.0*(np.random.uniform()-0.5) for i in range(len(t))]
counter = 0

for i in t:
  data.append(50.0*np.cos(2*np.pi*10.0*i)+40.0*np.cos(2*np.pi*20.0*i) + 30.0*np.cos(2*np.pi*30.0*i) + 20.0*np.cos(2*np.pi*40.0*i))
  counter+=1
print 'Size of signal: ',len(data)

coord = Vec_cpp()
zero  = Vec_cpp()
coord.extend(i for i in data)
zero.extend((0)*i for i in data)
naff = NAFF()

########## Modify default values #######

#naff.fmax=3
#naff.set_merit_function("minimize_RMS_time")
#naff.set_window_parameter(param, 'c')

tune_all=naff.get_f(coord,zero)
for i in tune_all:
  print i
