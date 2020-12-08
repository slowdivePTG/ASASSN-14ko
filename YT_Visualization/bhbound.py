import yt
import yt.units as u
from yt.fields.api import ValidateParameter
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--directory', '-d', dest='DIR',
                    help='data directory', type=str)
parser.add_argument(
    '--chk', help='checkpoint to plot, e.g., 0080', type=str, default='0080')
args = parser.parse_args()

os.chdir(args.DIR)
ds = yt.load("multitidal_hdf5_chk_{}".format(args.chk))

edata = np.loadtxt('extras.dat', dtype='float64')
G = u.gravitational_constant
gmpt = G * edata[6] * u.g

odata = np.loadtxt('pruned_sinks_evol.dat', dtype='float64')

if odata[0, 14] > odata[1, 14]:
    part_tag_pt = odata[0, 0]
else:
    part_tag_pt = odata[1, 0]
odata_pt = odata[np.where(odata[:, 0] == part_tag_pt)[0]]
odata_ob = odata[np.where(odata[:, 0] != part_tag_pt)[0]]

ptvec = odata_pt[:, 2:8]
obvec = odata_ob[:, 2:8]
boundvec = obvec
time = odata_pt[:, 1]
tindex = abs(time - ds.current_time.v).argmin()


def _bhbound(field, data):
    pos = np.zeros(data['x'].shape, dtype='float64')
    vel2 = pos.copy()
    for i, ax in enumerate(['x', 'y', 'z']):
        pos += (data[ax] - ptvec[tindex, i] * u.cm)**2.
    for i, ax in enumerate(['velx', 'vely', 'velz']):
        vel2 += (data[ax] - ptvec[tindex, i + 3] * u.cm / u.s)**2.
    pot = -gmpt / np.sqrt(pos) + 0.5 * vel2
    return (pot)


def _tfallback(field, data):
    pot = data['bhbound']
    return (gmpt * np.pi / np.sqrt(2) / (-pot)**(3 / 2))


def _selfbound(field, data):
    vel2 = np.zeros(data['x'].shape, dtype='float64')
    V = [data["velx"], data["vely"], data["velz"]]
    for i, ax in enumerate(['velx', 'vely', 'velz']):
        vel2 += (data[ax] - boundvec[tindex, i + 3] * u.cm / u.s)**2.

    pot = 0.5 * data['gpot'] + 0.5 * vel2
    return (pot)


yt.add_field(("gas", "bhbound"),
             function=_bhbound,
             units="erg/g",
             take_log=False,
             force_override=True,
             sampling_type="cell")

yt.add_field(("gas", "selfbound"),
             function=_selfbound,
             units="erg/g",
             take_log=False,
             force_override=True,
             sampling_type="cell")

yt.add_field(("gas", "tfallback"),
             function=_tfallback,
             units="day",
             take_log=False,
             force_override=True,
             sampling_type="cell")

da = ds.all_data().cut_region('obj["bhbound"] < obj["selfbound"]')
ds = yt.load("multitidal_hdf5_chk_{}".format(args.chk))
pp = yt.SlicePlot(ds, "z", "bhbound")
#pp.set_cmap(field="bhbound", cmap='turbo')
pp.zoom(3)
pp.save('tfallback_{}_{}.pdf'.format(args.DIR, args.chk))
