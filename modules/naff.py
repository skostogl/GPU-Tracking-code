from NAFF import *
from itertools import islice

def get_x(turns, particles):
    data_all = []
    for particle in (particles):
      data_x = Vec_cpp()
      data_x.extend ([( b.x [particle]) for b in turns])
    #data_x[0]=5
    #for b in turns:
    #  print b.x[0]
    #  quit()
      data_all.append(data_x)
    return data_all

def get_xp(turns, particles):
    data_all= []
    for particle in (particles):
      data_xp = Vec_cpp()
      data_xp.extend ([( b.xp [particle]) for b in turns])
      data_all.append(data_xp)
    return data_all

def get_y(turns, particles):
    data_all= []
    for particle in (particles):
      data_y = Vec_cpp()
      data_y.extend ([( b.y [particle]) for b in turns])
      data_all.append(data_y)
    return data_all

def get_yp(turns, particles):
    data_all= []
    for particle in (particles):
      data_yp = Vec_cpp()
      data_yp.extend ([( b.yp [particle]) for b in turns])
      data_all.append(data_yp)
    return data_all

def get_z(turns, particles):
    data_all= []
    for particle in (particles):
      data_z = Vec_cpp()
      data_z.extend ([( b.z [particle]) for b in turns])
      data_all.append(data_z)
    return data_all

def get_d(turns, particles):
    data_all= []
    for particle in (particles):
      data_d = Vec_cpp()
      data_d.extend ([( b.d [particle]) for b in turns])
      data_all.append(data_d)
    return data_all

def naff(data,data_p,second_half=False):
    tunes=[]
    numb_particle = 0
    for particle in (data):
      print "Particle ", numb_particle+1
      tune = NAFF_f1( data[numb_particle],data_p[numb_particle] )
      if (second_half == True):
          tune = 1-tune
      tunes.append(tune)
      print "Tune from NAFF: ",tune
      numb_particle = numb_particle + 1
    return tunes



