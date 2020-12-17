"""
local version of my dmde and dmdt making jYT script

experimented with a few different ways to smooth the data
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import simps, cumtrapz
from scipy.interpolate import splev, splrep, interp1d, splint
from scipy import stats
import math
import re
import os

# TODO write this code below that automatically formats things for the input directory
WRITE_TO_STARS_LIBRARY_INPUT = False

#exec(open("/Users/lawsmith/Dropbox/jYT/my_settings.py").read())
day = 86400.
yr = 3.154e7
G = 6.674e-8
M_sun = 1.988435e33
R_sun = 6.955e10
#mpl.rcParams['figure.figsize'] = (6,15)
#mpl.rcParams['font.size'] = 12
M_bh = 1e6*M_sun

RECENTER = True 

n0 = 40

s_linear = 0  #1e-7

NUM_X_PER_LOG_INTERVAL = 200  #10000
TIME_DAYS_EXTEND_MAX = 365*100 # todo maybe change this to 1000 yrs later

ds = [
##                                                to center dmde                                                dmdt                          slope
##                                     R*         minmax  il,   ir                name            s_dm_de       tmin  tmax             s      EXT
##d[0]                                 d[1]       d[2]    d[3]  d[4] d[5] d[6]    d[-6]            d[-5]         d[-4] d[-3]            d[-2]  d[-1]
#['m0.3_p1_b0.6_300k_0100.dat',       0.2814,   'min',  -0.2, 0.,                  'm0.3_t0.0',     50,           -99,  3.02e2,          6e-3,  0.85],
#['m0.3_p1_b0.7_300k_0100.dat',       0.2814,   'shoulders',  -1.0,0.,0.,1.0,      'm0.3_t0.0',     50,           -99,  2.5e2,           1e-4,  0.99],
#['m0.3_p1_b0.8_300k_0100.dat',       0.2814,   'min',  -0.2, 0.,                  'm0.3_t0.0',     100,          -99,  5.5e2,           2e-4,  0.99],
###['m0.3_p1_b0.9_300k_0100.dat',       0.2814,   'min',  -0.5, 0.,                  'm0.3_t0.0',     20,           1.1e1,  1.2e3,         5e-4,  0.6], 
#['m0.3_p1_b0.9_300k_lr15_0102.dat',   0.2814,   'min',  -1, 0.,                    'm0.3_t0.0',     20,           1.2e1,  1.2e3,         5e-4,  0.6],
#['m0.3_p1_b1.0_300k_0100.dat',       0.2814,   'min',  -0.5, 0.,                  'm0.3_t0.0',     10,           1.05e1,  2e2,          1e-4,  0.99],
#['m0.3_p1_b2.0_300k_0106.dat',        0.2814,   'min',  -0.5, 0.,                 'm0.3_t0.0',     50,           8e0,  5e2,             3e-3,  0.5],
##
#['m0.3_p11_b0.6_300k_0100.dat',      0.2989,   'min',  -0.2, 0.,                  'm0.3_t1.0',     100,           -99,  3.59e2,         2e-3,  0.99],
#['m0.3_p11_b0.7_300k_0100.dat',      0.2989,   'min',  -0.2, 0.,                  'm0.3_t1.0',     50,            -99,  2.8e2,          1e-4,  0.99],
#['m0.3_p11_b0.8_300k_0100.dat',      0.2989,   'shoulders',-1.0,-0.5,0.5,1.0,     'm0.3_t1.0',     10,            -99,  5.5e2,          1e-4,  0.99],
#['m0.3_p11_b0.9_300k_0100.dat',      0.2989,   'min',  -0.5, 0.,                  'm0.3_t1.0',     10,            1.25e1,  1.7e2,         3e-5,  0.7],
#['m0.3_p11_b1.0_300k_0100.dat',      0.2989,   'min',  -0.2, 0.,                  'm0.3_t1.0',     5,             1.2e1,  5e2,          5e-4,  0.3],
#['m0.3_p11_b2.1_300k_0103.dat',      0.2989,   'min',  -0.5, 0.,                  'm0.3_t1.0',     5,             1e1,  1e3,            2e-3,  0.8],
##
###choose['m0.5_p1_b0.6_300k_0060.dat',0.4452,  'shoulders', -1,0,0.5,1,           'm0.5_t0.0',     40,           2.95e1,  3e2,          1.5e-2, 0.99],
#['m0.5_p1_b0.6_300k_0100.dat',       0.4452,   'downup', -0.1, 0.1,0.3,           'm0.5_t0.0',     50,           3e1,  3.e2,            4.7e-2, 0.99],
#['m0.5_p1_b0.8_300k_0101.dat',        0.4452,   'min', -0.5, 0.5,                 'm0.5_t0.0',     100,           -99,  3.e2,            1e-3, 0.99],
#['m0.5_p1_b1.0_300k_0101.dat',         0.4452,   'min', -0.5, 0.5,                'm0.5_t0.0',     100,           -99,  1.3e3,            1e-3, 0.99],
#['m0.5_p1_b1.15_300k_0101.dat',        0.4452,   'min', -0.5, 0.5,                'm0.5_t0.0',     50,           1e1,  2e2,            1.e-4, 0.99],
#['m0.5_p1_b1.4_300k_0101.dat',       0.4452,   'shoulders', -1,0,0.5,1.5,         'm0.5_t0.0',     50,           9e0,  1.e3,            3e-3,   0.4],
#['m0.5_p1_b2.6_300k_0100.dat',        0.4452,   'min', -0.5, 0.5,                 'm0.5_t0.0',     50,           6e0,  8e2,             7e-3, 0.99],
##
#['m0.5_p38_b0.6_300k_0100.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     100,          3.1e1,  4.7e2,         2.5e-1,  0.99],
#['m0.5_p38_b0.8_300k_0100.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     90,           -99,  4e2,             5e-4,    0.99],
#['m0.5_p38_b1.0_300k_0100.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     50,           -99,  9.9e2,           5e-4,    0.99],
#['m0.5_p38_b1.2_300k_0100.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     5,            1.05e1,  8e2,          4e-4,    0.99],
#['m0.5_p38_b1.4_300k_0098.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     5,            1e1,  1.4e2,             1.1e-5,    0.956],
#['m0.5_p38_b2.8_300k_0106.dat',      0.4654,   'min',  -0.2, 0.,                  'm0.5_t1.0',     50,            7e0,  1.4e3,             1.25e-2,    0.6],
##
#['m0.7_p1_b0.75_300k_0060.dat',      0.6485,   'min',  -0.5,0.5,                  'm0.7_t0.0',     30,            2.5e1,  3.28e2,        8e-2,   0.8],
#['m0.7_p1_b1.0_300k_0060.dat',       0.6485,   'min',  -0.5,0.5,                  'm0.7_t0.0',     30,            1.35e1,  2.5e2,        2e-4,   0.8],
#['m0.7_p1_b1.25_300k_0100.dat',      0.6485,   'min',  -0.5,0.5,                  'm0.7_t0.0',     50,            1.2e1,  6.465e2,       4e-4,   0.8],
#['m0.7_p1_b1.5_300k_0100.dat',       0.6485,   'shoulders', -2,-1,0.5,1.5,        'm0.7_t0.0',     50,            8e0,  2.4e2,             5e-4,   0.99],
#['m0.7_p1_b3.4_300k_0100.dat',       0.6485,   'min', -0.5, 0.,                   'm0.7_t0.0',     100,             -99,  4.5e2,         1.2e-3, 0.7],
###choose['m0.7_p1_b3.4_300k_0075.dat',       0.6485,   'min', -0.5, 0.,           'm0.7_t0.0',     100,             5e0,  5e2,            2e-2, 0.75],
##
#['m0.7_p50_b0.8_300k_0060.dat',      0.6793,    'min',  -0.5, 0.,                 'm0.7_t1.0',     50,           2.4e1,  8e2,            8e-2,    0.6],
#['m0.7_p50_b1.0_300k_0060.dat',      0.6793,    'min',  -0.5, 0.,                 'm0.7_t1.0',     100,          1.55e1,   4.5e2,        3e-3,    0.6],
#['m0.7_p50_b1.15_300k_0100.dat',     0.6793,    'min',  -0.4, 0.,                 'm0.7_t1.0',     50,           1.7e1,   5.5e2,         3e-3,    0.6],
#['m0.7_p50_b1.3_300k_0100.dat',      0.6793,    'shoulders',  -1,-0.5,1,1.5,      'm0.7_t1.0',     100,            -99,  6e2,            7e-3,    0.6],
###choose['m0.7_p50_b1.3_300k_0100.dat',0.6793,  'max',  0,0.5,                    'm0.7_t1.0',     100,            -99,  1.5e3,          7e-3,    0.99],
#['m0.7_p50_b1.5_300k_0100.dat',      0.6793,    'max',  0.0,  0.4,                'm0.7_t1.0',     20,             -99,  1.5e3,          2e-3,    0.9],
#['m0.7_p50_b1.7_300k_0104.dat',      0.6793,    'shoulders',  -1.5,-1.,0.5,1,     'm0.7_t1.0',     10,             -99,  7e2,            7e-4,    0.8],
#['m0.7_p50_b3.6_300k_0075.dat',      0.6793,    'min',  0.5,  1.5,                'm0.7_t1.0',     50,             5e0,  4.5e2,            8.e-3,    0.4],
##
###choose['m1.0_p1_b1.0_300k_0040.dat',0.9012,   'shoulders',  -2.,0.,0.,2.,       'm1.0_t0.0',     15,           1.5e1,  3e2,            3e-3,     0.95],
#['m1.0_p1_b1.0_300k_0060.dat',       0.9012,    'min',  -0.5, 0.,                 'm1.0_t0.0',     50,           1.85e1,  1.49e2,        3e-2,     0.99],
#['m1.0_p1_b1.25_300k_0060.dat',      0.9012,    'min',  -0.5, 0.,                 'm1.0_t0.0',     50,           1.1e1,  2.5e2,          2e-3,     0.9],
#['m1.0_p1_b1.5_300k_0090.dat',       0.9012,    'min',  -0.4, 0.,                 'm1.0_t0.0',     50,            1e1,  6e2,             1e-3,     0.85],
#['m1.0_p1_b1.75_300k_0091.dat',      0.9012,    'min',  -0.4, 0.,                 'm1.0_t0.0',     20,           -99,  5e2,              5e-3,     0.5],
#['m1.0_p1_b2.0_300k_0081.dat',       0.9012,    'shoulders',  -2.,-1.,1.,1.5,     'm1.0_t0.0',     20,            6e0,  4e2,             2e-3,     0.6],
#['m1.0_p1_b4.2_300k_0080.dat',       0.9012,    'min',  0.5, 1.5,                 'm1.0_t0.0',     50,            -99,  1e3,             1.5e-2,     0.99],
##
###choose['m1.0_p10_b1.0_300k_0040.dat',1.0455,  'shoulders',   -1.5,-0.5,0.5,1.5,  'm1.0_t0.57',    50,            2e1,  2.5e2,           4.5e-2,    0.99],
#['m1.0_p10_b1.0_300k_0060.dat',      1.0455,    'shoulders',   -1.5,-0.5,0.5,1.5,  'm1.0_t0.57',    100,           2.6e1,  2.02e2,        3.3e-2,    0.99],
#['m1.0_p10_b1.5_300k_0060.dat',      1.0455,    'shoulders',   -1.5, -1., 1., 1.5, 'm1.0_t0.57',    50,            1.1e1,  6.05e2,        1.e-2,     0.99],
#['m1.0_p10_b2.0_300k_0075.dat',      1.0455,    'min',  -0.4, 0.,                  'm1.0_t0.57',    50,            1.1e1,  8e2,           5e-4,      0.99],
#['m1.0_p10_b2.5_300k_0061.dat',      1.0455,    'max',  0., 1.,                    'm1.0_t0.57',    10,            6e0,  4e2,             8e-4,      0.99],
#['m1.0_p10_b3.0_300k_0060.dat',      1.0455,    'shoulders',  -2.,-1.,1.,2.,       'm1.0_t0.57',    100,           3e0,  3.5e2,             4e-3,      0.99],
#['m1.0_p10_b3.5_300k_lr16_0062.dat',1.0455,     'shoulders',  -2,-1,1,2,           'm1.0_t0.57',    50,            3e0,  3e2,             2e-3,    0.6],
#['m1.0_p10_b4.9_300k_0062.dat',      1.0455,    'min',  -0.5, 0.5,                 'm1.0_t0.57',    100,           4e0,  3e2,             1e-3,      0.9],
##
#['m1.0_p16_b1.0_300k_lr16_0040.dat',  1.2872,  'min', -0.5, 0.5,                    'm1.0_t1.0',    100,           3e1, 2.57e2,           2e-1,     0.98],
###choose['m1.0_p16_b1.0_300k_lr16_0060.dat',  1.2872,  'min', -0.5, 0.5,            'm1.0_t1.0',    100,           3.6e1, 2.6e2,           3.5e-2,     0.9],
#['m1.0_p16_b1.5_300k_lr16_1110.dat',  1.2872,  'min', -0.5, 0.5,                    'm1.0_t1.0',    100,           1.8e1, 5e2,           1.e-1,     0.92],
#['m1.0_p16_b2.0_300k_lr16_0060.dat',1.2872,     'min', -0.5, 0.5,                   'm1.0_t1.0',    50,           1.65e1, 1.35e3,         6e-2,     0.9],
#['m1.0_p16_b3.0_300k_lr16_0062.dat', 1.2872,    'shoulders',   -2.5, -1.3, 1.2, 2.5,'m1.0_t1.0',    20,           9.5,  8e2,              2e-3,     0.9],
#['m1.0_p16_b4.0_300k_lr16_0052.dat', 1.2872,    'max',  0., 0.5,                    'm1.0_t1.0',    10,           7, 1e3,                 1.5e-3,   0.5],
#['m1.0_p16_b4.5_300k_lr16_0051.dat', 1.2872,    'max', 0., 1.,                      'm1.0_t1.0',    10,           7, 1e3,                 3e-3,     0.6],
#['m1.0_p16_b5.0_300k_lr16_0052.dat', 1.2872,    'max', 0., 1.,                      'm1.0_t1.0',    10,           7.5,   1e3,             1e-3,     0.99],
#['m1.0_p16_b6.0_300k_lr16_0054.dat', 1.2872,    'max', -1, 0,                       'm1.0_t1.0',    10,           8,   5e2,             5e-4,     0.6],
##
#['m1.5_p1_b1.0_300k_0060.dat',       1.6275,    'shoulders',  -2,0,0,2,             'm1.5_t0.0',    50,           3.5e1,    1.5e2,          2e-2,     0.99],
#['m1.5_p1_b1.5_300k_0060.dat',       1.6275,    'min',  -0.5, 0.5,                  'm1.5_t0.0',    50,           1.6e1,    2e2,          1.45e-1,     0.8],
###choose['m1.5_p1_b1.5_300k_0075.dat',       1.6275,    'min',  -0.5, 0.5,          'm1.5_t0.0',    50,           1.8e1,    1e3,          2.5e-1,     0.6],
#['m1.5_p1_b2.0_300k_2000_0076.dat',   1.6275,    'min',  -0.5, 0.5,                 'm1.5_t0.0',    50,           1.3e1,     1e3,          1.2e-4,     0.99],
#['m1.5_p1_b2.75_300k_0060.dat',      1.6275,    'max',  0,2,                        'm1.5_t0.0',    100,          6e0,      1.5e2,            1e-4,     0.99],
#['m1.5_p1_b6.7_300k_0055.dat',       1.6275,    'min',  -0.5, 0.5,                  'm1.5_t0.0',    50,           -99,      5e2,          4e-3,     0.8],
##
#['m1.5_p17_b2.0_300k_lr16_0060.dat',  2.0805,   'min',  -0.5, 0.5,                  'm1.5_t1.0',    100,           1.5e1,    8e2,         3e-1,    0.5],
#['m1.5_p17_b4.0_300k_lr16_0067.dat',   2.0805,   'max',  -0.5, 0.,                  'm1.5_t1.0',    10,            8e0,    8e2,           1e-3,    0.6],
#['m1.5_p17_b6.0_300k_lr16_0070.dat',  2.0805,   'max',  -0.5, 0.5,                  'm1.5_t1.0',    10,            9e0,    8e2,           3e-4,    0.5],
#['m1.5_p17_b8.6_300k_lr16_0055.dat',  2.0805,   'max',  -0.5, 0.5,                  'm1.5_t1.0',    10,            6e0,    2e3,           3e-3,    0.8],
##
#['m3.0_p16_b1.5_300k_lr16_0040.dat',    3.3192,     'min',  -0.5, 0.5,              'm3.0_t1.0',     30,            2.5e1,  4e2,           3.5e-2,    0.99],
#['m3.0_p16_b2.0_300k_lr16_0040.dat',    3.3192,     'min',  -0.5, 0.5,              'm3.0_t1.0',     30,            1.6e1,  5e2,           4e-2,    0.9],
###choose['m3.0_p16_b2.0_300k_lr16_0060.dat',3.3192, 'min',  -0.5, 0.5,              'm3.0_t1.0',     100,           1.6e1,  1e3,           2e-1,    0.6],
#['m3.0_p16_b3.0_300k_lr16_0065.dat',    3.3192,     'shoulders', -3,-1,1,3,         'm3.0_t1.0',     50,            1.2e1,  4e2,           1e-3,    0.5],
#['m3.0_p16_b4.0_300k_lr16_0053.dat',    3.3192,     'shoulders', -3,-1,1,3,         'm3.0_t1.0',     50,            1.1e1,  8.1e2,           1.1e-3,    0.99],
['m3.0_p16_b4.5_300k_lr14_0035.dat',    3.3192,     'max',  0., 0.25,               'm3.0_t1.0_lr14',     50,            1e1,  1e2,             5e-4,    0.8],
###choose['m3.0_p16_b4.5_300k_lr15_0052.dat', 3.3192,'max',  0., 0.6,                'm3.0_t1.0_lr15',     10,            1e1,  1.4e3,           1e-3,    0.9],
#['m3.0_p16_b4.5_300k_lr16_0052.dat',    3.3192,     'max',  0., 0.25,               'm3.0_t1.0',     50,            1e1,  5e2,             5e-4,    0.6],
#['m3.0_p16_b5.0_300k_lr16_0051.dat',    3.3192,     'shoulders', -3,-1,1,3,         'm3.0_t1.0',     50,             8e0,  5e2,            1e-3,    0.99],
#['m3.0_p16_b7.0_300k_lr16_0051.dat',    3.3192,     'max', 0,1,                     'm3.0_t1.0',     50,             9e0,  5e2,            1e-4,    0.5],
#['m3.0_p16_b10.8_300k_lr16_0058.dat',    3.3192,     'shoulders', -1,0,2,3,         'm3.0_t1.0',     50,             5.5e0,  5e2,            2.5e-3,    0.99],
##
#['m10.0_p3_b1.5_300k_0100.dat',       4.073,    'min',  -0.5, 0.,                   'm10.0_t0.1',   50,            1.1e1,  6e2,             4e-4,     0.99],
]

def find_bin_centers(x):
    c = []
    for i, xx in enumerate(x):
        if (i+1)==len(x):
            continue
        c.append((x[i]+x[i+1])/2.)
    return c

def round_up(n, decimals=0): 
    multiplier = 10 ** decimals 
    return math.ceil(n * multiplier) / multiplier

def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier) / multiplier

outputfile_integrated_dmde = open('/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/integrated_dmde_orig.dat', 'w')
outputfile_integrated_dmdt = open('/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/integrated_dmdt_Bspline_1e4yr.dat', 'w')

for d in ds:
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    #fig2, (ax2, ax3, ax4) = plt.subplots(3, 1)
    #fig3, ax3 = plt.subplots()
    #normalize to \Delta e, as in Ryu+2020b, Rees, Phinney, Stone?
    #delta_e = (G * M_star / R_star) * (M_bh / M_star)**(1./3)
    delta_e =  G * (float(d[0][1:4])*M_sun)**(2./3) * M_bh**(1./3) / (d[1]*R_sun)

    e, dm = np.loadtxt('/Users/lawsmith/Dropbox/1e6-grid/local_results/hists/'+d[0], skiprows=4)
    de = e[1]-e[0]
    # prune where dm=0, for b-spline fitting later
    sel = np.where((dm!=0.0) & np.isfinite(dm) & np.isfinite(e))
    e = e[sel]
    dm = dm[sel]
    dm_de = dm/de

    e_orig = e
    dm_de_orig = dm_de

    # real orig
    s_dm_de_orig = gaussian_filter(dm_de_orig, d[-5], mode='nearest')

    # shift dmde histogram to be centered at 0
    if RECENTER:
        if d[2] == 'min':
            e_orig = e_orig - e_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])][np.argmin(s_dm_de_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])])]
        elif d[2] == 'max':
            e_orig = e_orig - e_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])][np.argmax(s_dm_de_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])])]
        elif d[2] == 'shoulders':
            center = (e_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])][np.argmax(s_dm_de_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])])] \
                     + e_orig[(d[5]<e_orig/delta_e) & (e_orig/delta_e<d[6])][np.argmax(s_dm_de_orig[(d[5]<e_orig/delta_e) & (e_orig/delta_e<d[6])])]) / 2.
            e_orig = e_orig - center
        elif d[2] == 'shift':
            e_orig = e_orig + d[3]*delta_e
        elif d[2] == 'downup':
            center = (e_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])][np.argmin(s_dm_de_orig[(d[3]<e_orig/delta_e) & (e_orig/delta_e<d[4])])] \
                     + e_orig[(d[4]<e_orig/delta_e) & (e_orig/delta_e<d[5])][np.argmax(s_dm_de_orig[(d[4]<e_orig/delta_e) & (e_orig/delta_e<d[5])])]) / 2.
            e_orig = e_orig - center

    #ax.scatter(e/1e17, log_dm_de, c='r', rasterized=True, s=1, edgecolors='none')
    #ax.plot(e/1e17, slog_dm_de, c='g')
    #ax.scatter(e/delta_e, dm_de * delta_e/(float(d[0][1:4])*M_sun), c='r', rasterized=True, s=1, edgecolors='none')
    ax.scatter(e_orig/delta_e, dm_de_orig * delta_e/(float(d[0][1:4])*M_sun), c='k', rasterized=True, s=1, edgecolors='none')
    #ax.plot(e/delta_e, s_dm_de * delta_e/(float(d[0][1:4])*M_sun), c='k')
    ax.plot(e_orig/delta_e, s_dm_de_orig * delta_e/(float(d[0][1:4])*M_sun), c='r')
    #ax.axvline(ymin=-99,ymax=99,c='r')

    # calculate mdot
    e_bound_orig = e_orig[np.where(e_orig<0.)]
    t_orig = 2.*np.pi*G*M_bh/((2*np.abs(e_bound_orig))**(3./2))
    de_dt_orig = (1./3)*((2.*np.pi*G*M_bh)**(2./3))*t_orig**(-5./3)
    s_dm_de_bound_orig = s_dm_de_orig[np.where(e_orig<0.)]
    s_mdot_orig = s_dm_de_bound_orig*de_dt_orig

    dm_de_bound_orig = dm_de_orig[np.where(e_orig<0.)]
    mdot_orig = dm_de_bound_orig*de_dt_orig

    # integrate dmde or mdot as a new way to consistently calculate deltaM after shifting dmde
    #integrated_dmdt_orig = simps(mdot_orig, t_orig)
    integrated_dmde_orig = simps(dm_de_orig, e_orig)
    #outputfile_integrated_dmde.write(d[0] + '\t\t\t' + str(integrated_dmde_orig/M_sun) + '\t\t\t' + str(2*integrated_dmdt_orig/M_sun) + '\t\t\t' + str(2*integrated_dmdt_orig/integrated_dmde_orig) + '\n')
    outputfile_integrated_dmde.write(d[0] + '\t\t\t' + str(integrated_dmde_orig/M_sun) + '\n')

    # orig mdots
    ax2.scatter(t_orig/day, mdot_orig*yr/M_sun, c='k', rasterized=True, s=1, edgecolors='none')
    #ax3.scatter(t_orig/day, mdot_orig*yr/M_sun, c='k', rasterized=True, s=1, edgecolors='none')

    #changed both
    sel = np.where((t_orig/day > d[-4]) & (t_orig/day < d[-3]))
    t_trim = t_orig[sel]
    mdot_trim = mdot_orig[sel]

    # bin mdots
    # equally spaced bins
    # todo jamie come back to this
    log_range = np.log10(t_trim[-1]/day) - np.log10(t_trim[0]/day)
    my_bins = np.logspace(np.log10(t_trim[0]), np.log10(t_trim[-1]), num=int(round(log_range*n0, 0)))
    bin_centers = np.array(find_bin_centers(my_bins))
    y_binned = stats.binned_statistic(t_trim, mdot_trim, statistic='mean', bins=my_bins)

    # some y_binned[0] are nan because nothing in that bin
    sel = np.where(np.isfinite(y_binned[0]))
    t_trim = t_trim[sel]
    mdot_trim = mdot_trim[sel]
    y_binned_0 = y_binned[0][sel]
    bin_centers = bin_centers[sel]

    ax2.scatter(bin_centers/day, y_binned_0*yr/M_sun, c='b', s=10)

    spl = splrep(bin_centers/day, y_binned_0*yr/M_sun, s=s_linear)
    x2 = np.logspace(np.log10(min(bin_centers/day)), np.log10(max(bin_centers/day)), num=int(round(log_range*NUM_X_PER_LOG_INTERVAL, 0)))
    y2 = splev(x2, spl)
    ax2.plot(x2, y2, c='g')
    nts = (x2/y2)*splev(x2, spl, der=1)

    spllog = splrep(np.log10(bin_centers/day), np.log10(y_binned_0*yr/M_sun), s=d[-2])
    x2log = np.log10(x2)
    y2log = splev(x2log, spllog)
    #ax2.plot(10**x2log, 10**y2log, c='r')    
    ntslog = splev(x2log, spllog, der=1)

    x = 10**x2log
    y = 10**y2log

    ax2.axvline(max(x), c='k')

    # interpolate
    # jamie todo implement this x equally spaced in log with log_range
    #f = interp1d(x, y)
    log_range = np.log10(max(x))-np.log10(min(x))
    #xnew = np.logspace(np.log10(round_up(min(x), 8)), np.log10(round_down(max(x),8)), num=int(round(log_range,3)*1000))
    #ynew = f(xnew)
    #x = xnew
    #y = ynew

    # extend, n_inf, asymptotic power law index
    # todo temp
    slope_sel = np.where(x > d[-1]*max(x))
    ax2.axvline(x[slope_sel][0], c='k')
    slope = np.mean(ntslog[slope_sel])
    extend_log_range = np.log10(TIME_DAYS_EXTEND_MAX)-np.log10(max(x))
    extend_x = np.logspace(np.log10(max(x)), np.log10(TIME_DAYS_EXTEND_MAX), num=int(round(extend_log_range*NUM_X_PER_LOG_INTERVAL, 0)))[1:]
    b = np.log10(y[-1]) - slope*np.log10(x[-1])
    extend_y_log = slope*np.log10(extend_x) + b
    extend_y = 10**extend_y_log
    #ax2.plot(x, y, c='r')
    #ax2.plot(extend_x, extend_y, c='k')
    x = np.concatenate((x,extend_x))
    y = np.concatenate((y,extend_y))

    # n(t), power law index as a function of time
    #spl = splrep(np.log10(x), np.log10(y), k=3, s=0)
    #nt = splev(np.log10(x), spl, der=1)
    nt = np.concatenate((ntslog, np.full(len(extend_x), slope)))

    ascii.write([x,y],'/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/'+d[0],names=['t (day)','dm/dt (M_sun/yr)'],overwrite=True)
    if not os.path.exists('/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/' + d[-6]):
        os.makedirs('/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/' + d[-6])
    ascii.write([x,y],'/Users/lawsmith/Dropbox/1e6-grid/local_results/data-dmdts/' + d[-6] + '/' + d[0].split('_')[2][1:].ljust(5, '0') + '.dat', names=['t (day)','dm/dt (M_sun/yr)'],overwrite=True)

    ascii.write([[x[np.argmax(y)]],[max(y)]],'/Users/lawsmith/Dropbox/1e6-grid/local_results/data-tpeakmdotpeak/'+d[0],names=['t_peak (day)','{dm/dt}_peak (M_sun/yr)'],overwrite=True)
    ascii.write([[slope]],'/Users/lawsmith/Dropbox/1e6-grid/local_results/data-ninf/'+d[0],names=['n_inf'],overwrite=True)
    ascii.write([x,nt],'/Users/lawsmith/Dropbox/1e6-grid/local_results/data-nt/'+d[0],names=['t (day)','n(t)'],overwrite=True)

    ax2.plot(x, y, c='r')
    #ax2.scatter(x, y, c='r', s=0.1, edgecolors=None)

    #TODO 
    # extend x even further, to test integral of dmdt
    # todo I chose 1e4 years. probably OK.
    extend_log_range = np.log10(365*1e4)-np.log10(max(x))
    # todo maybe change this num depending on how long this takes to run
    extend_x = np.logspace(np.log10(max(x)), np.log10(365*1e4), num=100000)[1:]
    b = np.log10(y[-1]) - slope*np.log10(x[-1])
    extend_y_log = slope*np.log10(extend_x) + b
    extend_y = 10**extend_y_log
    xx = np.concatenate((x,extend_x))
    yy = np.concatenate((y,extend_y))

    # todo this is linear. ok?
    xx_cgs = xx*day
    yy_cgs = yy*M_sun/yr
    spl = splrep(xx_cgs, yy_cgs, k=3, s=0)
    integral_dmdt = splint(min(xx_cgs), max(xx_cgs), spl)
    outputfile_integrated_dmdt.write(d[0] + '\t\t\t' + str(2*integral_dmdt/M_sun) + '\n')


    # todo move / make general. fast overplot of GRR2013 5/3.
    #g13_53 = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/5-3/' + BETA + '.dat')
    # adjust to M~0.3 msun, R~0.28 rsun
    # M* = float(d[0][1:4]),     R* = d[1]
    #ax2.plot(g13_53[:,0] - np.log10(float(d[0][1:4])) + 1.5*np.log10(d[1]), g13_53[:,1] + 2.0*np.log10(float(d[0][1:4])) - 1.5*np.log10(d[1]), c='g', label='GRR13')

    ax.set_xlim(-5, 5)
    ax.set_xlabel(r'$e/\Delta e$')
    ax.set_ylabel(r'$dM/de\ [M_\star / \Delta e]$')
    ax.set_yscale('log')
    ax.set_title(re.sub('.dat', '', d[0]), fontsize=12)
    ax.grid(which='both')
    fig.tight_layout(pad=0.3)
    fig.savefig('/Users/lawsmith/Dropbox/1e6-grid/local_results/dmdes/dmde_'+re.sub('dat', '', d[0].replace(".","_"))+'_'+str(d[-5])+'.pdf')

    ax2.set_xlim(1., 5e4)
    #ax2.set_ylim(1e-7, 1.2e1)
    ax2.set_ylim(1e-7, 1e2)
    #ax2.set_xlim(1e1, 1e3)
    #ax2.set_ylim(1e-4, 1e1)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    #ax2.set_title(re.sub('.dat', '', d[0]), fontsize=12)
    ax2.set_ylabel(r'$\dot M\ {\rm [M_\odot/yr]}$')
    ax2.set_xlabel(r'$t\ \mathrm{[day]}$')
    ax2.grid(which='both')
    #fig2.tight_layout(pad=0.3)
    #fig2.savefig('/Users/lawsmith/Dropbox/1e6-grid/local_results/dmdts/dmdt_'+re.sub('dat', '', d[0].replace(".","_"))+'_'+str(d[-2])+'.pdf')

    ax3=ax2.twinx()
    #ax3.plot(x2, nts, c='g')
    ax3.axhline(-2.2, c='gray', ls='-')
    ax3.plot(x, nt, c='k')
    #ax3.set_ylim(-5, 15)
    ax3.set_ylim(-4, 4)
    ax3.axhline(-5/3, c='gray', ls='--')
    ax3.axhline(-9/4, c='gray', ls='--')

    ax4=ax2.twinx()
    ax4.plot(x2, cumtrapz(y2, x2*day/yr, initial=0)*2/(integrated_dmde_orig/M_sun), c='g')
    #ax4.plot(10**x2log, cumtrapz(10**y2log, (10**x2log)*day/yr, initial=0)*2/(integrated_dmde_orig/M_sun), c='r')
    ax4.plot(x, cumtrapz(y, x*day/yr, initial=0)*2/(integrated_dmde_orig/M_sun), c='k')
    ax4.set_ylim(0,1.1)
    ax4.axhline(1, c='k')
    ax4.set_yticks((0,1))

    fig2.tight_layout(pad=0.3)
    fig2.savefig('/Users/lawsmith/Dropbox/1e6-grid/local_results/dmdts_combo/'+re.sub('dat', '', d[0].replace(".","_"))+'_'+str(d[-2])+'.pdf')

    """
    ax3.set_xlim(0, 100)
    ax2.set_ylim(1e-5, 1e1)
    ax3.set_ylim(bottom=1e-4)
    ax3.set_yscale('log')
    ax3.set_title(re.sub('.dat', '', d[0]), fontsize=12)
    ax3.set_ylabel(r'$\dot M\ {\rm [M_\odot/yr]}$')
    ax3.set_xlabel(r'$t\ \mathrm{[day]}$')
    ax3.grid(which='both')
    #ax3.grid(which='minor', lw=0.1, ls=':')
    fig3.tight_layout(pad=0.3)
    fig3.savefig('/Users/lawsmith/Dropbox/1e6-grid/local_results/dmdts_linear/dmdt_linear_'+re.sub('dat', '', d[0].replace(".","_"))+'_'+str(d[-2])+'.pdf')
    """
    plt.close()

outputfile_integrated_dmde.close()
outputfile_integrated_dmdt.close()