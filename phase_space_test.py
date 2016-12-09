import pickle
from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

lattice = Lattice()
if ( True ):
#if ( False ):
  lattice.read_twiss_table("LHC/lhc_no_bb.twi")
  lattice.optimise()
  lattice.compile()
  lattice.write_ptx("LHC/lhc_no_bb")
else:
  lattice.read_ptx("LHC/lhc_no_bb")

lattice.n_turns = 5000
lattice.norm_emit_x = 2e-6
lattice.norm_emit_y = 2e-6
lattice.collect_tbt_data = 1 # every 1 turn

#b,grid = cmp_grid (lattice.sigma_x(), lattice.sigma_x()*8, lattice.sigma_y(), lattice.sigma_y()*8,1,lattice)
#n_particles=b.size()
n_particles=20
b=HostBunch(n_particles)
for i in range (n_particles):
    b.y[i]=0
    b.yp[i]=0
    b.x[i]=(i+1)*lattice.sigma_x()
    b.xp[i]=0
    #b.d[i]=100e-5


#lattice.track(b)
#for i in range (n_particles):
#  filename = '/home/skostogl/cuTrack/dat_files/particles/particle_%d.dat'%i
#  tbt = [ (b.x[i], b.xp[i], b.y[i], b.yp[i]) for b in lattice.turns ]
#  with open(filename,'w') as outfile:
#    for t in tbt:
#      outfile.write("{} {} {} {}\n".format(t[0], t[1], t[2], t[3]))
#
#
#fig,ax=plt.subplots()
#color=iter(cm.rainbow (np.linspace(0,1,n_particles)))
#for i in range (n_particles):
#  c=next(color)
#  for j in (lattice.turns):
#    plt.plot(j.x[i]*1e3, j.xp[i]*1e3, c=c,marker='o',ms=3,markeredgewidth=0.0)
#
plt.xlabel(r'$x [mm]$', fontsize=20)
plt.ylabel(r'$x_p [mm]$', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
annotation_string=r'$\delta=0$'
at = AnchoredText(annotation_string,prop=dict(size=18), frameon=True,loc=1)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)
plt.tight_layout()
plt.show()
