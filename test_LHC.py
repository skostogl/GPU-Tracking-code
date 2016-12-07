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
  lattice.read_twiss_table("LHC/lhc_no_bb.twi")
  #lattice.read_twiss_table("LHC/lhc_Qp0.twi")
  #lattice.read_twiss_table("LHC/lhc_Qp0.twi")
  lattice.compile()
  lattice.write_ptx("LHC/lhc_no_bb")
  #lattice.write_ptx("LHC/lhc_Qp0")
else:
  lattice.read_ptx("LHC/lhc_no_bb")
  #lattice.read_ptx("LHC/lhc_Qp0")

lattice.n_turns = 1000
lattice.norm_emit_x = 2e-6
lattice.norm_emit_y = 2e-6
lattice.collect_tbt_data = 1 # every 1 turn

#print lattice.sigma_x(), lattice.sigma_y()
#b,grid = cmp_grid (lattice.sigma_x(), lattice.sigma_x()*4, lattice.sigma_y(), lattice.sigma_y()*4,0.3,lattice)
#n_particles=b.size()
n_particles=10
#b=HostBunch(n_particles)
b=lattice.make_matched_bunch(n_particles)
#for i in range (n_particles):
#    b.y[i]=0
#    b.yp[i]=0
#    b.d[i]=10e-5
lattice.track(b)

############################ FMA  ##################################
tunes_x = naff(lattice.turns[0:1000], vec_HostBunch.x, vec_HostBunch.xp, second_half=False)
tunes_y = naff(lattice.turns[0:1000], vec_HostBunch.y, vec_HostBunch.yp, second_half=False)
fig,ax=create_plot(tunes_x,tunes_y,0)
plt.show()

#tunes_x1, tunes_y1, tunes_x2, tunes_y2, tune_diffusion = FMA(lattice.turns[0:500], lattice.turns[500:1000])
#fig,ax=create_plot(tunes_x2,tunes_y2, 0, tune_diffusion, colorbar=True, resonance_diagram=False, order=4)
#plt.show()

#particle_x=[lattice.turns[0].x[i] for i in range (n_particles)]
#particle_y=[lattice.turns[0].y[i] for i in range (n_particles)]
#fig,ax = create_plot (particle_x,particle_y,0,tune_diffusion,colorbar=True,resonance_diagram=False,order=4)
#plt.show()



#for i in range (n_particles):
#  filename = '/home/skostogl/cuTrack/dat_files/test/particle_%d_pos_delta.dat'%i
#  tbt = [ (b.x[i], b.xp[i], b.y[i], b.yp[i]) for b in lattice.turns ]
#  with open(filename,'w') as outfile:
#    for t in tbt:
#      outfile.write("{} {} {} {}\n".format(t[0], t[1], t[2], t[3]))

filename = '/home/skostogl/cuTrack/dat_files/footprint_grid.dat'
tbt = [ (tunes_x[i],tunes_y[i]) for i in range (n_particles)]
with open(filename,'w') as outfile:
    for t in tbt:
      outfile.write("{} {} \n".format(t[0], t[1]))
##useful
#particles_x=[(lattice.turns[0].x[i]/lattice.sigma_x()) for i in range (n_particles)]
#particles_y=[(lattice.turns[0].y[i]/lattice.sigma_y()) for i in range (n_particles)]
#create_plot(particles_x,particles_y,grid)

#filename = '/home/skostogl/cuTrack/dat_files/initial_footprint_grid.dat'
#tbt = [ (particles_x[i],particles_y[i]) for i in range (n_particles)]
#with open(filename,'w') as outfile:
#    for t in tbt:
#      outfile.write("{} {} \n".format(t[0], t[1]))
#plt.show()


filename = '/home/skostogl/cuTrack/dat_files/grid_footprint_grid.dat'
tbt = [ (grid[i][0],grid[i][1]) for i in range (len(grid))]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {}\n".format(t[0], t[1]))
