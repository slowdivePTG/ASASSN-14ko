import yt
import yt.units as u
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from data_binning import *
plt.rcParams.update({
    "text.usetex": True})

G = u.gravitational_constant
q0 = 1e6
beta = 1
M = q0 * u.Msun
rT = q0**(1 / 3) * u.Rsun
deltaE = (G * M / rT**2 * u.Rsun).in_cgs()


class TDE:
    def __init__(self,
                 q=q0,
                 Ecc=1,
                 beta=1,
                 file='1.000.dat',
                 ms=1 * u.Msun,
                 rs=1 * u.Rsun,
                 DIR='../../STARS_library/retrieval/m1.0_t0.445/'):
        G = u.gravitational_constant
        q0 = 1e6
        M0 = q0 * u.Msun

        with open(DIR + file) as f:
            lines = (line for line in f if not line.startswith('"'))
            dmdt_t0 = np.loadtxt(lines, skiprows=1)

        t0 = dmdt_t0[:, 0]
        dmdt0 = dmdt_t0[:, 1]
        dedt0 = (2 / 3 * (G**2 * M**2 * np.pi**2 / 2)**(1 / 3) *
                 (t0 * u.day)**(-5 / 3)).in_cgs()
        dmde0 = (dmdt0 * (u.Msun / u.yr) / dedt0).in_cgs()

        rT = (q0 / (ms / u.Msun))**(1 / 3) * rs
        deltaE = (G * M0 / rT**2 * rs).in_cgs()
        E = -((G**2 * M0**2 * np.pi**2 / 2)**(1 / 3) *
              (t0 * u.day)**(-2 / 3)).in_cgs()
        E_dE = E / deltaE
        dmde_dE = dmde0 * deltaE

        self.M = q * u.Msun
        self.dM_Ms = (dmde0[1:] * (E[1:] - E[:-1])).sum() * 2 / ms
        self.Ecc = Ecc
        self.beta = beta
        self.rT = (self.M / ms)**(1 / 3) * rs
        self.deltaE = (G * self.M / self.rT**2 * rs).in_cgs()
        self.E = self.E_Ecc(E_dE * self.deltaE)
        self.T = ((-self.E**3 * 2)**(-1 / 2) * np.pi * G *
                  self.M).in_units('day')
        self.Tfallback = self.T[0]

        self.dmde = dmde_dE / self.deltaE
        dedt = 2 / 3 * ((G**2 * self.M**2 * np.pi**2 / 2)**(1 / 3) *
                        self.T**(-5 / 3)).in_cgs()
        self.dmdt = (dedt * self.dmde).in_units('Msun/yr')

        # while np.argmax(self.dmdt) + 1 == len(self.dmdt):
        while self.dmdt[-1] > self.dmdt[-2]:
            self.dmdt = self.dmdt[:-1]
            self.T = self.T[:-1]
        self.Tpeak = self.T[np.argmax(self.dmdt)].in_units('day')
        self.Tpeak0 = self.Tpeak - self.Tfallback

        if self.Eorb == 0:
            self.Period = np.inf
        else:
            self.Period = 2 * np.pi * G * self.M / (-2 * self.Eorb)**(3 / 2)
            self.Period = self.Period.in_units('day')

    def E_Ecc(self, E):
        if self.Ecc == 1:
            self.Eorb = 0
        else:
            a = self.rT / self.beta / (1 - self.Ecc)
            self.Eorb = -G * self.M / 2 / a
        return E + self.Eorb


def mesa_param(Dir):
    if Dir[4] == '_':
        Ms = float(Dir[1:4])
        age = float(Dir[6:])
    else:
        Ms = float(Dir[1:5])
        age = float(Dir[7:])
    if Ms == 0.1:
        if age == 0:
            rs = 0.1214
        else:
            rs = 0.1215
        rhoc_rho = 5.5
    elif Ms == 0.3:
        if age == 0:
            rs = 0.2814
        else:
            rs = 0.2989
        rhoc_rho = 5.8
    elif Ms == 0.5:
        if age == 0:
            rs = 0.4452
            rhoc_rho = 11
        else:
            rs = 0.4564
            rhoc_rho = 12
    elif Ms == 0.7:
        if age == 0:
            rs = 0.6485
            rhoc_rho = 23
        else:
            rs = 0.6793
            rhoc_rho = 36
    elif Ms == 1.0:
        if age == 0:
            rs = 0.9012
            rhoc_rho = 42
        elif age == 1:
            rs = 1.2872
            rhoc_rho = 756
        else:
            rs = 1.0455
            rhoc_rho = 138
    elif Ms == 1.5:
        if age == 0:
            rs = 1.6275
            rhoc_rho = 128
        else:
            rs = 2.0805
            rhoc_rho = 1697
    elif Ms == 3.0:
        if age == 0:
            rs = 1.8896
            rhoc_rho = 73
        else:
            rs = 3.3192
            rhoc_rho = 1198
    elif Ms == 10:
        if age == 0:
            rs = 3.6870
            rhoc_rho = 38
        else:
            rs = 8.4132
            rhoc_rho = 1292
    return Ms * u.Msun, rs * u.Rsun, rhoc_rho
