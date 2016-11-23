import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *

############################ LATTICE ##################################

lattice = Lattice()
lattice.BETX=0.4
lattice.BETY=0.4
lattice.energy=6500
lattice.norm_emit_x=2e-6
lattice.norm_emit_y=2e-6
lattice.add("BeamBeam(1.15e8,%d, 11e-6)" % lattice.rel_gamma())
#lattice.add("Multipole(3,1e8,0)")
lattice.optimise()
lattice.compile()
lattice.n_turns = 1
lattice.collect_tbt_data = 1 

############################ BUNCH  ##################################

b,grid = cmp_grid(lattice.sigma_x(),lattice.sigma_x()*6,lattice.sigma_y(),lattice.sigma_y()*6,0.5)
n_particles=b.size()
for turn in range(1000):
  r2 = b.x[0] *b.x[0]+b.y[0]*b.y[0]
  lattice.track(b)
  del lattice.turns[-1]
  for i in range(n_particles):
    x  = b.x[i]
    xp = b.xp[i]
    y  = b.y[i]
    yp = b.yp[i]
    phi_x = 0.31*6.28
    b.x[i]  =  x*np.cos(phi_x) + lattice.BETX * xp*np.sin(phi_x)
    b.xp[i] = -x*np.sin(phi_x) / lattice.BETX + xp*np.cos(phi_x)
    phi_y = 0.32*6.28
    b.y[i]  =  y*np.cos(phi_y) + lattice.BETY * yp*np.sin(phi_y)
    b.yp[i] = -y*np.sin(phi_y) / lattice.BETY + yp*np.cos(phi_y)
filename = 'tbt.dat'
tbt = [ (b.x[0], b.xp[0], b.y[0], b.yp[0], b.z[0] , b.d[0]) for b in lattice.turns ]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {} {} {} {} {}\n".format(t[0], t[1], t[2], t[3], t[4], t[5]))

############################ NAFF  ##################################

tunes_x = naff(lattice.turns[1:1000], vec_HostBunch.x, vec_HostBunch.xp)
tunes_y = naff(lattice.turns[1:1000], vec_HostBunch.y, vec_HostBunch.yp)
fig,ax=create_plot(tunes_x,tunes_y, grid,resonance_diagram=False)
plt.show()

############################ FMA  ##################################

#tunes_x1, tunes_y1, tunes_x2, tunes_y2, tune_diffusion = FMA(lattice.turns[1:500], lattice.turns[500:1000])
#fig,ax=create_plot(tunes_x2,tunes_y2, grid, tune_diffusion, colorbar=True, resonance_diagram=True, order=4)
#plt.show()


