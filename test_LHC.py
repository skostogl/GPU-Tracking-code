import pickle
from cuTrack import *

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
lattice.collect_tbt_data = 1 # every 1 turn

b = HostBunch(1)

lattice.track(b)

print lattice.turns[1000].x[0]

# save turn by turn x, xp data for particle 0
filename = 'tbt.dat'
tbt = [ (b.x[0], b.xp[0], b.y[0], b.yp[0]) for b in lattice.turns ]
with open(filename,'w') as outfile:
  for t in tbt:
    outfile.write("{} {} {} {}\n".format(t[0], t[1], t[2], t[3]))

############# FFT TUNE ##################

import matplotlib.pyplot as plt
import numpy as np

t = np.arange(lattice.n_turns)
sp = np.fft.fft( [lattice.turns[i].x[0] for i in range(lattice.n_turns)] )
freq = np.fft.fftfreq(t.shape[-1])
plt.plot(freq, sp.real, freq, sp.imag)

plt.show()

