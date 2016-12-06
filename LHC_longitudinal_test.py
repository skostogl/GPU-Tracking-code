import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *

lattice = Lattice()
if ( True ):
#if ( False ):
  #lattice.read_twiss_table("LHC/lhc_no_bb.twi")
  lattice.read_twiss_table("LHC/lhc_no_bb.twi")
  lattice.optimise()
  lattice.compile()
  #lattice.write_ptx("LHC/lhc_no_bb")
  lattice.write_ptx("LHC/lhc_no_bb")
else:
  #lattice.read_ptx("LHC/lhc_no_bb")
  lattice.read_ptx("LHC/lhc_no_bb")

lattice.n_turns = 1000
lattice.norm_emit_x = 2e-6
#lattice.bunch_length=1e-3
lattice.norm_emit_y = 2e-6

n_particles=50
b=lattice.make_matched_bunch(n_particles)

lattice.collect_tbt_data = 1 # every 1 turn

lattice.bunch_energy_spread = 1e-4

lattice.track(b)

############################ DELTA-Z  ##################################

for i in range (n_particles):
  print 'Particle: ',i
  plt.plot(lattice.turns[:].z(i), lattice.turns[:].d(i))
plt.show()

############################ NAFF  ##################################

tunes = naff(lattice.turns[0:1000], vec_HostBunch.z, vec_HostBunch.d)
