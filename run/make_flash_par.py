import mesa_reader as mr
import numpy as np
from yt.units import *


def make_flash_par(num):
    h = mr.MesaData('profile{}.data'.format(num))
    R = (h.photosphere_r * Rsun).in_cgs()
    M = (h.star_mass * Msun).in_cgs()
    tdyn = np.sqrt(R**3 / gravitational_constant / M).in_cgs()
    Trelax = tdyn * 5
    Tsim = tdyn * 100
    Mbh = 1e6 * Msun
    P = 110 * day
    a = ((gravitational_constant * Mbh / 4 / np.pi**2) * P**2)**(1 / 3)
    rT = (R * (Mbh / M)**(1 / 3)).in_cgs()
    beta = 0.6
    r_peri = rT / beta
    e = 1 - r_peri / a
    print(a.in_units('AU'), r_peri.in_units('AU'))

    with open('./flash.par', 'r') as f:
        lines = f.readlines()
        lines[3] = '# M = {:.4f} msun\n'.format(float(h.star_mass))
        lines[4] = '# R = {:.4f} rsun\n'.format(float(h.photosphere_r))
        lines[5] = '# tdyn = {:.4e} s\n'.format(float(tdyn))


        lines[28] = 'xmax = {:.6e}  # 1000 rstar\n'.format(float(R) * 1000)
        lines[29] = 'ymax = {:.6e}  # 1000 rstar\n'.format(float(R) * 1000)
        lines[30] = 'zmax = {:.6e}  # 500 rstar\n'.format(float(R) * 500)
        lines[31] = 'sim_xCenter = {:.6e}  # 500 rstar\n'.format(float(R) * 500)
        lines[32] = 'sim_yCenter = {:.6e}  # 500 rstar\n'.format(float(R) * 500)
        lines[33] = 'sim_zCenter = {:.6e}  # 250 rstar\n'.format(float(R) * 250)

        lines[21] = 'checkpointFileIntervalTime = {:.6e}   # 1 tdyn\n'.format(float(tdyn))
        lines[23] = 'plotFileIntervalTime       = {:.6e}   # 1 tdyn\n'.format(float(tdyn))
        lines[48] = 'tmax = {:.6e}   # 100 tdyn\n'.format(float(tdyn) * 100)
        lines[78] = 'sim_tRelax = {:.6e}   # 5 tdyn\n'.format(float(tdyn) * 5)

        lines[62] = 'sim_parentPeri = {:.8e}\n'.format(float(r_peri))
        lines[67] = 'sim_orbEcc    = {:.8e}\n'.format(float(e))
        lines[167] = 'sink_softening_radius = {:.8e} # ~1/2 rp\n'.format(float(r_peri)/2)

    with open('./flash.par', 'w') as f:
        for l in lines:
            f.write(l)


import sys
make_flash_par(sys.argv[1])