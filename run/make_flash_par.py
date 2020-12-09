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
    Mbh = 7e7 * Msun
    P = 110 * day
    #a = ((gravitational_constant * Mbh / 4 / np.pi**2) * P**2)**(1 / 3)
    rT = (R * (Mbh / M)**(1 / 3)).in_cgs()
    beta = 0.6
    r_peri = rT / beta
    #e = 1 - r_peri / a
    #print(a.in_units('AU'), r_peri.in_units('AU'))

    with open('./flash.par', 'r') as f:
        lines = f.readlines()
        lines[4] = '# M = {:.2f} msun\n'.format(float(h.star_mass))
        lines[5] = '# R = {:.2f} rsun\n'.format(float(h.photosphere_r))
        lines[6] = '# tdyn = {:.4e} s\n'.format(float(tdyn))

        lines[17] = 'checkpointFileIntervalTime = {:.2f}   # 1 tdyn\n'.format(
            float(tdyn))
        lines[19] = 'plotFileIntervalTime       = {:.2f}   # 100 tdyn\n'.format(
            float(tdyn * 100))

        lines[24] = 'xmax = {:.6e}  # 1000 rstar\n'.format(float(R) * 1000)
        lines[25] = 'ymax = {:.6e}  # 1000 rstar\n'.format(float(R) * 1000)
        lines[26] = 'zmax = {:.6e}  # 500 rstar\n'.format(float(R) * 500)
        lines[28] = 'sim_xCenter = {:.6e}  # 500 rstar\n'.format(
            float(R) * 500)
        lines[29] = 'sim_yCenter = {:.6e}  # 500 rstar\n'.format(
            float(R) * 500)
        lines[30] = 'sim_zCenter = {:.6e}  # 250 rstar\n'.format(
            float(R) * 250)

        lines[45] = 'tmax = {:.2f}   # 100 tdyn\n'.format(float(tdyn) * 100)
        lines[80] = 'sim_tRelax = {:.2f}   # 5 tdyn\n'.format(float(tdyn) * 5)

        lines[58] = 'sim_ptMass = {:.4e}\n'.format(float(Mbh.in_cgs()))
        lines[66] = 'sim_periBeta  = {:.1f}\n'.format(beta)
        lines[68] = 'sim_period    = {:.6e}\n'.format(float(P.in_units('s')))
        lines[170] = 'sink_softening_radius = {:.8e} # ~1/2 rp\n'.format(
            float(r_peri) / 2)

    with open('./flash.par', 'w') as f:
        for l in lines:
            f.write(l)


import sys
make_flash_par(sys.argv[1])
