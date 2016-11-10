import sys
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *
from modules.naff import *
from modules.grid import *


############################ LATTICE ##################################
lattice = Lattice()
lattice.read_twiss_table("small_ring/lattice_RF.twi")
lattice.optimise()
lattice.compile()
lattice.n_turns = 200
lattice.collect_tbt_data = 1 # every 1 turn
lattice.norm_emit_x=1e-5
lattice.norm_emit_y=1e-5
lattice.bunch_energy_spread= 1e-100
lattice.bunch_length=5.5e-3
############################ BUNCH  ##################################
n_particles=10
b=lattice.make_matched_bunch(n_particles)

for i in range (n_particles):
    b.x[i]=1e-3
    b.xp[i]=1e-3
    b.y[i]=1e-3
    b.yp[i]=1e-3
#n_particles=b.size()
#b, grid = cmp_grid(lattice.sigma_x,lattice.sigma_x*10,lattice.sigma_y,lattice.sigma_y*10,0.0002)
#create_plot(b.x,b.y,grid)   # grid optional parameter
#create_plot(b.x,b.y)
lattice.track(b)

filename = 'tbt.dat'
tbt = [ (b.x[0], b.xp[0], b.y[0], b.yp[0], b.z[0] , b.d[0]) for b in lattice.turns ]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {} {} {} {} {}\n".format(t[0], t[1], t[2], t[3], t[4], t[5]))

############################ FFT TUNE  ##################################
data_x  = get_x  (lattice.turns[:], particles = range(0,n_particles)) 
data_xp = get_xp (lattice.turns[:], particles = range(0,n_particles))
data_y  = get_y  (lattice.turns[:], particles = range(0,n_particles))
data_yp = get_yp (lattice.turns[:], particles = range(0,n_particles))
data_z  = get_z  (lattice.turns[:], particles = range(0,n_particles))
data_d  = get_d  (lattice.turns[:], particles = range(0,n_particles))

tunes_x = naff(data_x, data_xp, flag=1)
tunes_y = naff(data_y, data_yp, flag=1)
tunes_z = naff(data_z, data_d)


create_plot (data_z,data_d,title="z-d",xlabel="z",ylabel="d")
#create_plot (data_x,data_xp)
#create_plot (tunes_x,tunes_y)
#create_plot (tunes_x,tunes_y,grid)
