from tracker import *
import math
import time 
import itertools
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from matplotlib import collections as mc

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

def create_plot(data_x, data_y, index_list=0,title="",xlabel="",ylabel="" ):
  segment=[]
  fig,ax=plt.subplots()
  if index_list==0:    
    for i in range (len(data_x)):  
      ax.plot(data_x[i],data_y[i],marker='o',ms=10,alpha=1,color='b')
  else:
      
    for element in index_list:
      segment.append( [ ( data_x[element[0]],data_y[element[0]] ) ,( data_x[element[1]],data_y[element[1]] ) ])
          
    lc = mc.LineCollection(segment, colors='k', linewidths=2)
    ax.add_collection(lc)
    ax.autoscale()
  date   = time.strftime( "%Y-%m-%d-%H:%M:%S" )
  plt. xlabel (xlabel)
  plt. ylabel (ylabel)
  plt. title  (title)
  plt.show()
  #fig.savefig( 'saved_plots/'+date+'_'+title+'.pdf' )  

def cmp_grid(sigma_x_init, sigma_x_final, sigma_y_init, sigma_y_final, step): 
  width =  sigma_x_final - sigma_x_init + step
  height = sigma_y_final - sigma_y_init + step
  n_particles_x = int (math.ceil(width/step))
  n_particles_y = int (math.ceil(height/step))
  n_particles   = int (n_particles_x * n_particles_y)
  b = HostBunch(n_particles)
  particle = 0
  i = sigma_x_init
  while (i < sigma_x_final + step):
    j = sigma_y_init
    while (j < sigma_y_final + step ):
      b.x [particle] = i 
      b.y [particle] = j
      particle = particle + 1
      j = j + step
    i = i + step
  segment_indexes=cmp_index_segment(n_particles_x, n_particles_y)
  return b, segment_indexes


