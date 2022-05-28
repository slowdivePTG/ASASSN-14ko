import mesa_reader as mr
import numpy as np
from yt.units import *
import os
import argparse

parser = argparse.ArgumentParser(
    description='Rewrite the flash.par after renewing the stellar model and orbital parameters.')
parser.add_argument('--profile',
                    help='Profile number')
parser.add_argument('--period',
                    help='Orbital period (day)', default=-1, type=int)
parser.add_argument('--beta',
                    help='Penatration parameter', default=1.0, type=float)
args = parser.parse_args()


def make_flash_par(num, p, beta):
    h = mr.MesaData('profile{}.data'.format(num))
    R = (h.photosphere_r * Rsun).in_cgs()
    M = (h.star_mass * Msun).in_cgs()
    tdyn = np.sqrt(R**3 / gravitational_constant / M).in_cgs()
    Trelax = tdyn * 5
    Tsim = tdyn * 100
    Mbh = 7e7 * Msun
    #a = ((gravitational_constant * Mbh / 4 / np.pi**2) * P**2)**(1 / 3)
    rT = (R * (Mbh / M)**(1 / 3)).in_cgs()
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
        if p != -1:
            lines[68] = 'sim_period    = {:.6e}\n'.format(p * 86400)
        else:
            lines[68] = 'sim_orbEcc    = 0.99999\n'
        lines[170] = 'sink_softening_radius = {:.8e} # ~1/2 rp\n'.format(
            float(r_peri) / 2)

    with open('./flash.par', 'w') as f:
        for l in lines:
            f.write(l)


make_flash_par(num=args.profile, p=args.period, beta=args.beta)
os.system('rm sm.dat')
os.system(
    'python ~/jYT/mesa_wrapper.py -i profile{}.data -o sm.dat -e h1 he3 he4 li7 be9 c12 n14 o16 f19 ne20 na23 mg24 mg26 al27 si28 p31 s32'.format(args.profile))
