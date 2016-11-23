import sys
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *

############################ LATTICE ##################################

lattice = Lattice()
lattice.read_twiss_table("small_ring/lattice_RF.twi")
lattice.optimise()
lattice.compile()
lattice.n_turns =1000
lattice.collect_tbt_data = 1 # every 1 turn
lattice.norm_emit_x=1e-5
lattice.norm_emit_y=1e-5
lattice.bunch_energy_spread=1e-4
#lattice.bunch_length=1e-3

############################ BUNCH  ##################################
n_particles=50
b=lattice.make_matched_bunch(n_particles)

lattice.track(b)

filename = 'tbt.dat'
tbt = [ (b.x[0], b.xp[0], b.y[0], b.yp[0], b.z[0] , b.d[0]) for b in lattice.turns ]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {} {} {} {} {}\n".format(t[0], t[1], t[2], t[3], t[4], t[5]))

############################ DELTA-Z  ##################################

for i in range (n_particles):
  print 'Particle: ',i
  plt.plot(lattice.turns[:].z(i), lattice.turns[:].d(i))
plt.show()

############################ NAFF  ##################################

tunes = naff(lattice.turns[0:1000], n_particles, vec_HostBunch.z, vec_HostBunch.d)
