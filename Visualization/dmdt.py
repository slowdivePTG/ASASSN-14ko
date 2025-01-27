import numpy as np
import yt
import yt.units as u
import matplotlib.pyplot as plt
import glob

from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import splev, splrep, interp1d, splint
###   splrep(x, y[, w, xb, xe, k, task, s, t, …])   ###
# Find the B-spline representation of a 1-D curve.

###   splev(x, tck[, der, ext])   ###
# Evaluate a B-spline or its derivatives.

###   splint(a, b, tck[, full_output])   ###
# Evaluate the definite integral of a B-spline between two given points.

from scipy import stats
import os
import argparse
import sys


def bins(arr, n, type='mean'):
    def mode(arr):
        hist, bin_edges = np.histogram(arr, bins=20)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        return bin_centers[np.argmax(hist)]
    bins = len(arr) // n
    if type == 'mean':
        arr2 = np.array([arr[i * n:i * n + n].mean() for i in range(bins)])
    elif type == 'median':
        arr2 = np.array([np.median(arr[i * n:i * n + n])
                         for i in range(bins)])
    elif type == 'mode':
        arr2 = np.array([10**mode(np.log10(arr[i * n:i * n + n]))
                         for i in range(bins)])
    return arr2


'''
def bins_above(x, y, n0, n1=3, thresh=5e-3):
    x = bins(x, n0)
    y = bins(y, n0)
    flag = False
    cut = []
    for i in range(len(x)):
        if (not flag) and y[i] > thresh:
            flag = True
            cut.append(i)
        if (flag) and y[i] < thresh:
            flag = False
            cut.append(i)
    if len(cut) % 2 != 0:
        cut.append(-1)
    elif len(cut) == 0:
        return x, y
    cut.append(-1)
    X, Y = x[:cut[0]], y[:cut[0]]
    for i in range(len(cut) // 2):
        X = np.append(X, bins(x[cut[2 * i]:cut[2 * i + 1]], n1))
        X = np.append(X, x[cut[2 * i + 1]:cut[2 * i + 2]])
        Y = np.append(Y, bins(y[cut[2 * i]:cut[2 * i + 1]], n1))
        Y = np.append(Y, y[cut[2 * i + 1]:cut[2 * i + 2]])
    return X, Y
'''


def find_bin_centers(x):
    return (x[1:] + x[:-1]) / 2


def smooth(x_trim, y_trim, n=100, NUM_X_PER_INTERVAL=1000, log=True):
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
    spl = splrep(bin_centers, y_binned_0, s=1e-2)
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
        spl = splrep(bin_centers, y_binned_0)
        y2 = splev(x2, spl)
        nts = splev(x2, spl, der=1)
        x = x2
        y = y2
    return (x, y)


class dmdt:
    def __init__(self, DIR, chk, Period=-1, Ecc=None,
                 bin=10000, D=188 * u.Mpc, eta=0.06, ax=None,
                 Mh=7e7 * u.Msun, M=1 * u.Msun, R=1 * u.Rsun):
        os.chdir('/Users/chang/Desktop/Santa Cruz/TDE_plot')
        os.chdir(DIR)

        if Period < 0:
            if Ecc == None:
                filename = 'b{}_dm_de_chk_{}.npz'.format(bin, chk)
            else:
                filename = 'b{}_dm_de_chk_{}_*.npz'.format(bin, chk)
                filename = glob.glob(filename)[0]
        else:
            filename = 'b{}_dm_de_chk_{}_{}.npz'.format(bin, chk, Period)

        npzfile = np.load(filename)
        output = npzfile['x']
        output_bd = npzfile['y']
        output_mdot = npzfile['z']
        self.DIR = DIR
        self.beta = float(self.DIR[10:13])
        rT = (Mh / M)**(1 / 3) * R
        rp = rT / self.beta
        if Period < 0:
            self.Period = np.inf * u.day
        else:
            self.Period = Period * u.day
        if Ecc == None:
            self.Ecc = float(1 - (4 * np.pi**2 * rp**3 / self.Period **
                                  2 / u.gravitational_constant / Mh)**(1 / 3))
        else:
            self.Ecc = Ecc
        self.e = output[:, 0]
        self.dm_de, self.s_dm_de = output[:, 1], output[:, 2]

        self.t = output_mdot[:, 0]
        self.mdot, self.mdot_orig = output_mdot[:, 1], output_mdot[:, 2]

        self.F = (self.mdot * u.Msun / u.yr * eta * u.c**2 / np.pi / 4 /
                  D**2).in_units('erg/s/cm**2')
        self.F_orig = (self.mdot_orig * u.Msun / u.yr * eta * u.c**2 / np.pi /
                       4 / D**2).in_units('erg/s/cm**2')
        if Period < 0:
            self.label = 'e = {:.1f}'.format(self.Ecc)
        else:
            self.label = 'e = {:.4f} (P = {})'.format(self.Ecc, self.Period)

    def dm_de_e(self, ax, Bin=100, n=60,
                color='grey', linestyle='-',
                Mh=7e7 * u.Msun,
                M=1 * u.Msun, R=1 * u.Rsun, beta=1,
                shift=False,
                search_range=0.05):
        # dm_de v.s. e
        try:
            if ax == None:
                f, ax = plt.subplots(figsize=(8, 6))
        except:
            pass
        s_e, s_dm_de = bins(self.e, n=50, type='mean'), bins(
            self.dm_de, n=50, type='mean')
        deltae = (u.gravitational_constant * Mh / R *
                  (M / Mh)**(2 / 3)).in_cgs()  # * beta**2
        dm_de = (self.dm_de * u.g**2 / u.erg * deltae).in_units('Msun')

        if shift:
            search = np.argwhere(abs(s_e / deltae) < search_range)
            arg = np.argmin(s_dm_de[search])
            e_shift = (s_e / deltae)[search][arg].v
        else:
            e_shift = 0
        dEde = (self.e / deltae).v - e_shift
        dEde_bin = np.linspace(-2, 2, Bin)
        width = dEde_bin[1] - dEde_bin[0]
        dm_de_binned = np.array([
            dm_de[abs(dEde - dEde_bin[i]) < width / 2].mean() for i in range(Bin)])

        # binned
        ax.scatter(dEde_bin,
                   dm_de_binned,
                   # linestyle=linestyle,
                   s=10,
                   color=color,
                   alpha=0.8)

        # smooth
        '''arg = np.argwhere(dm_de_binned > 1e-2).flatten()
                                dEde_spline, s_dm_de_spline = smooth(
                                    dEde_bin[arg[0] - 1:arg[-1] +
                                             1], dm_de_binned[arg[0] - 1:arg[-1] + 1],
                                    n=25,
                                    NUM_X_PER_INTERVAL=1000,
                                    log=False)
                                ax.plot(dEde_spline,
                                        s_dm_de_spline,
                                        linewidth=3,
                                        # linestyle=linestyle,
                                        color=color,
                                        label=self.label)'''

        # original
        ax.plot(dEde, dm_de, alpha=0.2, color=color)
        ax.set_yscale('log')
        ax.set_xlabel(r'$E/\delta E$', fontsize=25)
        ax.set_ylabel(
            r'$\mathrm dM$/$\mathrm dE$ ($M_\odot/\delta E$)', fontsize=25)
        ax.tick_params(labelsize=20)

    def Mdot_t(self, ax, Flux=False, norm=False, normfactor=1, Bin=1000, N=30,
               Mh=7e7 * u.Msun, M=1 * u.Msun, R=1 * u.Rsun):
        # Mdot v.s. t
        try:
            if ax == None:
                f, ax = plt.subplots(figsize=(8, 6))
        except:
            pass
        s_t, s_F, s_mdot = bins(self.t, Bin), bins(
            self.F, Bin), bins(self.mdot, Bin)
        s_t, s_F = smooth(self.t, self.F, n=N,
                          NUM_X_PER_INTERVAL=Bin, log=True)
        s_t, s_mdot = smooth(self.t, self.mdot, n=N,
                             NUM_X_PER_INTERVAL=Bin, log=True)
        peak_index = np.argmax(s_mdot)
        arg_efold = np.argwhere(s_mdot >= s_mdot.max() / np.e).flatten()
        self.efolding_rise = s_t[peak_index] - s_t[arg_efold[0]]
        self.efolding_fall = s_t[arg_efold[-1]] - s_t[peak_index]
        if norm:
            tindex = np.argmax(s_F)
            ax.plot((s_t - s_t[tindex]) * normfactor,
                    s_F / np.max(s_F),
                    label=self.label)
            ax.set_ylabel(r'Normalized Flux', fontsize=25)
        elif Flux:
            # ax.plot(self.t, self.F_orig, alpha=0.3)
            ax.plot(s_t,
                    s_F,
                    label=self.label)
            ax.set_ylabel(r'$F$ (erg/s/cm$^2$)', fontsize=25)
        else:
            # ax.plot(self.t, self.mdot_orig, alpha=0.3)
            ax.plot(s_t,
                    s_mdot,
                    label=self.label)
            ax.set_ylabel(r'$\dot M$ (Msun/yr)', fontsize=25)

        ax.set_xlabel('Time (day)', fontsize=25)
        ax.tick_params(labelsize=20)
        plt.tight_layout()
