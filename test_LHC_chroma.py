import pickle
from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *
import timeit
import matplotlib.pyplot as plt
import numpy as np

lattice = Lattice()
if ( True ):
#if ( False ):
  lattice.read_twiss_table("LHC/lhc_Qp0.twi")
  lattice.optimise()
  lattice.compile()
  lattice.write_ptx("LHC/lhc_Qp0")
else:
    lattice.read_ptx("LHC/lhc_Qp0")

lattice.n_turns = 10000
lattice.norm_emit_x = 2e-6
lattice.norm_emit_y = 2e-6
lattice.collect_tbt_data = 1 # every 1 turn

print lattice.sigma_x(), lattice.sigma_y()
b,grid = cmp_grid (lattice.sigma_x(), lattice.sigma_x()*6, lattice.sigma_y(), lattice.sigma_y()*6,0.03)
n_particles=b.size()

#for i in range (n_particles):
#    b.y[i]=0
#    b.yp[i]=0
#    b.d[i]=65e-5

lattice.track(b)

############################ FMA  ##################################
tunes_x = naff(lattice.turns[0:10000], vec_HostBunch.x, vec_HostBunch.xp, second_half=False)
tunes_y = naff(lattice.turns[0:10000], vec_HostBunch.y, vec_HostBunch.yp, second_half=False)
fig,ax=create_plot(tunes_x,tunes_y)
plt.show()

tunes_x1, tunes_y1, tunes_x2, tunes_y2, tune_diffusion = FMA(lattice.turns[0:5000], lattice.turns[5000:10000])
fig,ax=create_plot(tunes_x2,tunes_y2, 0, tune_diffusion, colorbar=False, resonance_diagram=False, order=4)
plt.show()
filename = '/home/skostogl/cuTrack/dat_files/FMA_chroma0.dat'
tbt = [(tunes_x2[i], tunes_y2[i],tune_diffusion[i]) for i in range (n_particles)]
with open(filename,'w') as outfile:
    for t in tbt:
      outfile.write("{} {} {}\n".format(t[0], t[1], t[2]))

particle_x=[(lattice.turns[0].x[i]/lattice.sigma_x()) for i in range (n_particles)]
particle_y=[(lattice.turns[0].y[i]/lattice.sigma_y()) for i in range (n_particles)]
fig,ax = create_plot (particle_x,particle_y,0,tune_diffusion,colorbar=True,resonance_diagram=False,order=4)
plt.show()

