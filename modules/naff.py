from NAFF import *
import numpy as np
#from multiprocessing import Pool
#from multiprocessing.dummy import Pool as ThreadPool
def naff(data, coord, coord_prime=0, second_half =False):
  tunes=[]
  n_particles = len(data[0].x)
  for i in range (n_particles):
    #tune=NAFF_f1(coord(data,i),coord_prime(data,i))
    tune=NAFF.get_f1(coord(data,i),coord_prime(data,i))
    if (second_half == True):
      tune = 1-tune
    print 'NAFF for particle %i and for %i turns : %f' %(i+1,len(data), tune)
    tunes.append(tune)
  return tunes

#def naff(data, coord, coord_prime=0, second_half =False):
#  n_particles = len(data[0].x)
#  pool=Pool(processes=8)
#  #tunes=[pool.apply_async(NAFF_f1,(coord(data,i),coord_prime(data,i))) for i in range (n_particles)]
#  tunes1, tunes2, tunes3=zip(*pool.map(NAFF_f1, ))
#  return tunes
