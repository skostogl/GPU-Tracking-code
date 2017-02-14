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

lattice.n_turns = 1000
lattice.norm_emit_x = 2e-6
lattice.norm_emit_y = 2e-6
lattice.collect_tbt_data = 1 # every 1 turn

n_particles=10
b=HostBunch(n_particles)

d=np.arange(0,35e-5,35e-5/n_particles)
for i in range (n_particles):
    b.x[i]=lattice.DX*d[i]
    b.y[i]=lattice.DY*d[i]
    b.xp[i]=lattice.DPX*d[i]
    b.yp[i]=lattice.DPY*d[i]
    b.d[i]=d[i]
print lattice.DX,lattice.DPX,lattice.DY,lattice.DPY

lattice.track(b)

tunes = naff(lattice.turns[0:1000], vec_HostBunch.z, vec_HostBunch.d)
fig,ax=plt.subplots()
color=iter(cm.rainbow (np.linspace(0,1,n_particles)))
for i in range (n_particles):
  c=next(color)
  for j in (lattice.turns):
    plt.plot(j.z[i]*1e3, j.d[i]*1e3, c=c,marker='o',ms=3,markeredgewidth=0.0)


plt.xlabel(r'$z [10^{-3}]$', fontsize=20)
plt.ylabel(r'$\delta [10^{-3}]$', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title("Longitudinal phase space", fontsize=20)
plt.tight_layout()
plt.show()

