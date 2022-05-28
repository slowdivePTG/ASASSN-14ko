import sys
import numpy as np
from yt.units import cm, second
filename = sys.argv[1]
ds = yt.load(filename)
da = ds.all_data()

#for i in ds.derived_field_list:
    #print(i)

print('CFL limit')
cmax = np.max([da['sound_speed'].in_cgs(), da['velocity_x'].in_cgs()], axis=0)*cm/second
temp = (da['dx']/cmax).in_cgs()
temp.sort()
print(temp[:5])

print('Particle crossing limit')
print('Particle velocity x', da['particle_velx'].in_units('km/s'))
print(0.5*(da['dx'].min()/da['particle_velx']).in_cgs())
print('Particle velocity y', da['particle_vely'].in_units('km/s'))
print(0.5*(da['dy'].min()/da['particle_vely']).in_cgs())
print('Particle velocity z', da['particle_velz'].in_units('km/s'))
print(0.5*(da['dz'].min()/da['particle_velz']).in_cgs())
