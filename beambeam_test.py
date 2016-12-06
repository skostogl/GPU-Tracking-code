import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *

############################ LATTICE ##################################

lattice = Lattice()
lattice.BETX = 1
lattice.BETY = 1
lattice.energy = 2
lattice.norm_emit_x = 1e-5
lattice.norm_emit_y = 1e-5
lattice.mass = 0.000510998928
print lattice.rel_gamma(), lattice.energy
lattice.add("BeamBeam(1.15e8,%d, 5e-5)" % lattice.energy)
lattice.optimise()
lattice.compile()
lattice.n_turns = 1000
lattice.collect_tbt_data = 1

############################ BUNCH  ##################################

initial_conditions = np.arange(-5*lattice.sigma_x(), 5*lattice.sigma_x(), step=lattice.sigma_x()/10)
n_particles = len(initial_conditions)
b = HostBunch(n_particles)
for i in range (n_particles):
    b.xp[i] = 0
    b.x[i] = initial_conditions[i]
    b.yp[i] = 0
    b.y[i] = initial_conditions[i]
print lattice.sigma_x()
lattice.track(b)

filename = 'tbt.dat'
tbt = [ (b.x[0], b.xp[0], b.y[0], b.yp[0], b.z[0] , b.d[0]) for b in lattice.turns ]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {} {} {} {} {}\n".format(t[0], t[1], t[2], t[3], t[4], t[5]))

############################ FMA  ##################################

fig, ax =create_plot (lattice.turns[-1].x, lattice.turns[-1].xp)
plt.show()
