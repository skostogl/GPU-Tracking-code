import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *

# LATTICE 
lattice = Lattice()
lattice.add("Drift(1)")
lattice.add("Quad(0.3)")
lattice.add("Drift(1)")
lattice.add("Quad(-0.3)")

lattice.compile()
lattice.n_turns =1000
lattice.collect_tbt_data = 1 # every 1 turn

# BUNCH 
bunch = HostBunch(10)
bunch.particle[0].x = 1

# RUN
lattice.track(bunch)

# PLOT
plt.plot(lattice.turns.x(0), lattice.turns.xp(0), '.')
plt.show()

