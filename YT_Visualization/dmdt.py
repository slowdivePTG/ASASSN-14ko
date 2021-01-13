import numpy as np
import yt
import yt.units as u
import matplotlib.pyplot as plt

from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import splev, splrep, interp1d, splint
from scipy import stats
import os
import argparse
import sys


def bins(arr, n):
    arr2 = [arr[i:i + n + 1].mean() for i in range(len(arr) - n)]
    return arr2


def find_bin_centers(x):
    return (x[1:] + x[:-1]) / 2


def smooth(x_trim, y_trim, n=100, NUM_X_PER_INTERVAL=1000, log=False):
    if log:
        log_range = np.log10(x_trim[-1] / x_trim[0])
        my_bins = np.logspace(np.log10(x_trim[0]), np.log10(
            x_trim[-1]), num=int(round(n * log_range, 0)))
    else:
        linear_range = x_trim[-1] - x_trim[0]
        my_bins = np.linspace(x_trim[0], x_trim[-1], num=int(round(n, 0)))
    bin_centers = np.array(find_bin_centers(my_bins))
    y_binned = stats.binned_statistic(
        x_trim, y_trim, statistic='mean', bins=my_bins)

    # some y_binned[0] are nan because nothing in that bin
    sel = np.where(np.isfinite(y_binned[0]))
    x_trim = x_trim[sel]
    y_trim = y_trim[sel]
    y_binned_0 = y_binned[0][sel]
    bin_centers = bin_centers[sel]
    spl = splrep(bin_centers, y_binned_0, s=0)
    if log:
        x2 = np.logspace(np.log10(min(bin_centers)),
                         np.log10(max(bin_centers)),
                         num=int(round(NUM_X_PER_INTERVAL, 0)))
        y2 = splev(x2, spl)
        spllog = splrep(np.log10(bin_centers), np.log10(y_binned_0), s=None)
        x2log = np.log10(x2)
        y2log = splev(x2log, spllog)
        ntslog = splev(x2log, spllog, der=1)
        x = 10**x2log
        y = 10**y2log
    else:
        x2 = np.linspace(min(bin_centers),
                         max(bin_centers),
                         num=int(round(NUM_X_PER_INTERVAL, 0)))
        y2 = splev(x2, spl)
        spl = splrep(bin_centers, y_binned_0, s=None)
        y2 = splev(x2, spl)
        nts = splev(x2, spl, der=1)
        x = x2
        y = y2
    return (x, y)


class dmdt:
    def __init__(self, DIR, chk, Period, D=188 * u.Mpc, eta=0.06, ax=None):
        os.chdir('/Users/chang/Desktop/Santa Cruz/TDE_plot')
        os.chdir(DIR)

        if Period < 0:
            filename = 'b10000_dm_de_chk_{}.npz'.format(chk)
        else:
            filename = 'b10000_dm_de_chk_{}_{}.npz'.format(chk, Period)

        npzfile = np.load(filename)
        output = npzfile['x']
        output_bd = npzfile['y']
        output_mdot = npzfile['z']
        self.DIR, self.Period = DIR, Period
        self.e = output[:, 0]
        self.dm_de, self.s_dm_de = output[:, 1], output[:, 2]

        self.t = output_mdot[:, 0]
        self.mdot, self.mdot_orig = output_mdot[:, 1], output_mdot[:, 2]

        self.F = (self.mdot * u.Msun / u.yr * eta * u.c**2 / np.pi / 4 /
                  D**2).in_units('erg/s/cm**2')
        self.F_orig = (self.mdot_orig * u.Msun / u.yr * eta * u.c**2 / np.pi /
                       4 / D**2).in_units('erg/s/cm**2')
        if Period < 0:
            self.label = self.DIR
        else:
            self.label = '{}, P = {}'.format(self.DIR, self.Period)

    def dm_de_e(self, ax, Bin=100, color='grey', Mh=7e7 * u.Msun, M=1 * u.Msun, R=1 * u.Rsun, beta=1):
        # dm_de v.s. e
        try:
            if ax == None:
                f, ax = plt.subplots(figsize=(8, 6))
        except:
            pass
        s_e, s_dm_de = bins(self.e, n=Bin), bins(self.dm_de, n=Bin)
        # s_e, s_dm_de = smooth(self.e, self.dm_de, n=Bin)
        deltae = (u.gravitational_constant * Mh / R *
                  (M / Mh)**(2 / 3)).in_cgs() * beta**2
        ax.plot(s_e / deltae,
                s_dm_de,
                # s=1,
                color=color,
                label='{}, P = {}'.format(self.DIR, self.Period))
        ax.plot(self.e / deltae, self.dm_de, alpha=0.2, color=color)
        ax.set_yscale('log')
        ax.set_xlabel(r'$\epsilon/\Delta\epsilon$', fontsize=20)
        ax.set_ylabel(r'd$M$/d$\epsilon$ (g$\cdot$s$^2$/cm$^2$)', fontsize=20)

    def Mdot_t(self, ax, Flux=False, norm=False, normfactor=1, Bin=100, N=30):
        # Mdot v.s. t
        try:
            if ax == None:
                f, ax = plt.subplots(figsize=(8, 6))
        except:
            pass
        s_t, s_F = smooth(self.t, self.F_orig, n=N,
                          NUM_X_PER_INTERVAL=Bin, log=True)
        s_t, s_mdot = smooth(self.t, self.mdot_orig, n=N,
                             NUM_X_PER_INTERVAL=Bin, log=True)
        if norm:
            tindex = np.argmax(s_F)
            ax.plot((s_t - s_t[tindex]) * normfactor,
                    s_F / np.max(s_F),
                    label='{}, P = {}'.format(self.DIR, self.Period))
            ax.set_ylabel(r'Normalized Flux', fontsize=0)
        elif Flux:
            ax.plot(self.t, self.F_orig, alpha=0.3)
            ax.plot(self.t,
                    self.F,
                    label='{}, P = {}'.format(self.DIR, self.Period))
            ax.set_ylabel(r'$F$ (erg/s/cm$^2$)', fontsize=20)
        else:
            ax.plot(self.t, self.mdot_orig, alpha=0.3)
            ax.plot(s_t,
                    s_mdot,
                    label='{}, P = {}'.format(self.DIR, self.Period))
            ax.set_ylabel(r'$\dot M$ (Msun/yr)', fontsize=20)

        ax.set_xlabel('Time (day)', fontsize=20)
        plt.tight_layout()
