from tracker import *
import math
import itertools
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from tune_resonances import *

def cmp_index_segment(n_particles_x, n_particles_y):
  n_particles=int (n_particles_x * n_particles_y)
  indexes = []
  index_list  = []
  index_list_unq = []
  for index in range (n_particles):
    indexes=[ item for item in (index-1, index-n_particles_y, index+n_particles_y, index+1 ) if item >=0 and item < n_particles ]
    if  ( (index+1)%(n_particles_y) == 0):
      del indexes[-1:]
    elif ( index%(n_particles_y) == 0):
      del indexes[0]
    index_list.extend( list(product([index], indexes)) )
  index_list = [ sorted(item) for item in (index_list) ]
  index_list.sort()
  index_list_unq = list(index_list for index_list ,_ in itertools.groupby(index_list))
  return index_list_unq

def create_plot(data_x, data_y,index_list=0, diff_tunes=0, colorbar=False,resonance_diagram=False,order=0 ):
  segment=[]
  fig,ax=plt.subplots()
  if (index_list == 0):    
    if (colorbar == False):  
      for i in range (len(data_x)):  
        plt.plot(data_x[i], data_y[i], marker='o',ms=3, alpha=1, color='b')
  else:
    for element in index_list:
      segment.append( [ ( data_x[element[0]],data_y[element[0]] ) , ( data_x[element[1]],data_y[element[1]] ) ])
    lc = mc.LineCollection(segment, colors='navy', linewidths=1)
    ax.add_collection(lc)
    plt.plot()
  if (colorbar == True and diff_tunes != 0):
    x = data_x
    y = data_y
    z = diff_tunes
    plt.scatter(x, y, edgecolors='none', s=15, c=np.log10(z))
    cb =plt.colorbar()
    cb.set_label(r'$log\sqrt{\Delta Q^2_x + \Delta Q^2_y}$', fontsize=25)
    cb.ax.tick_params(labelsize=20)
    plt.plot()
  if (resonance_diagram == True):
    if (order==0):
      print "Maximum order of resonance_diagram not defined!!"
    else:
      qx = np.sum(data_x)/len(data_x)
      qy = np.sum(data_y)/len(data_y)
      make_resonance_diagram(order, [0,1], [0,1])
      plt.plot()
  plt.axis('equal')
  return fig, ax


def cmp_grid(sigma_x_init, sigma_x_final, sigma_y_init, sigma_y_final, step,lattice): 
  n_particles_x =  int(round((sigma_x_final-sigma_x_init)/(step*lattice.sigma_x())))
  print n_particles_x
  if (sigma_x_init+n_particles_x*step*lattice.sigma_x()+step*lattice.sigma_x()<sigma_x_final):
      n_particles_x=n_particles_x+1
  #n_particles_x =  int(round((sigma_x_final)/(step*sigma_x_init)))
  #n_particles_x = int (round(sigma_x_final/(sigma_x_init*step)))
  n_particles_y=n_particles_x
  n_particles   = int (n_particles_x * n_particles_y)
  b = HostBunch(n_particles)
  particle = 0
  #for i in range (1,n_particles_x+1):
  #  for j in range (1,n_particles_y+1):
  for i in range (n_particles_x):
    for j in range (n_particles_y):
      #b.x [particle] =  i*step*sigma_x_init
      #b.y [particle] =  j*step*sigma_y_init
      b.x [particle] = sigma_x_init+ i*step*lattice.sigma_x()
      b.y [particle] = sigma_y_init+ j*step*lattice.sigma_y()
      particle = particle + 1
  #particles_x=[(b.x[i]/lattice.sigma_x() ) for i in range (particle)]
  #particles_y=[(b.y[i]/lattice.sigma_y() ) for i in range (particle)]
  #create_plot(particles_x,particles_y)
  segment_indexes=cmp_index_segment(n_particles_x, n_particles_y)
  return b, segment_indexes


