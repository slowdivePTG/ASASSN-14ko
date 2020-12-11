import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt.units as u
import yt
import os

part_tag = 0
time = 1
posx, posy, posz = 2, 3, 4
velx, vely, velz = 5, 6, 7
accelx, accely, accelz = 8, 9, 10
anglx, angly, anglz = 11, 12, 13
mass, mdot, ptime = 14, 15, 16


class SinkEvol:
    def __init__(self, DIR, runs, force, clean):
        self.user = DIR
        self.runs = runs
        self.Tdyn = 0

        import paramiko
        dirs = os.listdir('../../TDE_plot/')

        filename = '../../TDE_plot/sinks_evol.dat'
        if DIR != 'ptgcliu':
            runs = '{}_{}'.format(DIR[:3], self.runs)
        self.output = '../../TDE_plot/{}_pruned_sinks_evol.dat'.format(
            runs)
        if clean:
            os.system('rm ../../TDE_plot/{}*'.format(runs))

        exist = False
        for i in dirs:
            temp = '{}_pruned_sinks_evol.dat'.format(runs)
            if temp in i:
                exist = True
                break

        if not exist or force:
            hostname = "fend02.hpc.ku.dk"
            port = 22
            username = "ptgcliu"

            client = paramiko.SSHClient()
            key = paramiko.RSAKey.from_private_key_file(
                '/Users/chang/.ssh/id_rsa', password='I miss Shiyu')
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(hostname, port, username, pkey=key, compress=True)
            sftp_client = client.open_sftp()
            try:
                sftp_client.get(
                    '/storage/dark/{}/{}/extras.dat'.format(self.user, self.runs), '../../TDE_plot/{}_extras.dat'.format(runs))
            except FileNotFoundError:
                print('{} not found!'.format(
                    '/storage/dark/{}/{}/extras.dat'.format(self.user, self.runs)))
                os.system('rm ../../TDE_plot/{}_extras.dat'.format(runs))

            file_cluster = '/storage/dark/{}/{}/sinks_evol.dat'.format(
                DIR, self.runs)
            # print(file_cluster)
            try:
                sftp_client.get(file_cluster, filename)
            except FileNotFoundError:
                print('{} not found!'.format(file_cluster))
                return
            df = pd.read_table(filename, delim_whitespace=True)
            # search for headers printed again in a restart
            args = df[df['[01]time'] == '[01]time'].index

            # stop when no header's left
            while len(args) != 0:
                # the last restart
                restart_time = df['[01]time'][args[-1] + 1]
                df_before = df[df.index < args[-1]]
                # check if the restart time occured before
                arg_temp = df_before[df_before['[01]time']
                                     == restart_time].index

                # remove all data between two same restart time
                if len(arg_temp) > 0:
                    remove = range(arg_temp[0], args[-1] + 1)
                # otherwise remove only the header
                else:
                    remove = [args[-1]]
                df = df.drop(remove)
                args = df[df['[01]time'] == '[01]time'].index
            df.to_csv(self.output,
                      sep=' ',
                      index=False,
                      header=False)
        df = np.loadtxt(self.output)
        part_tag_bh = df[0, part_tag] if df[0, mass
                                            ] > df[1, mass] else df[1, part_tag]

        self.bh = df[df[:, part_tag] == part_tag_bh, :]
        self.star = df[df[:, part_tag] != part_tag_bh, :]
        self.T = self.star[:, time]
        os.system('rm ../../TDE_plot/sinks_evol.dat')
        self.runs = runs

        extras = np.loadtxt('../../TDE_plot/{}_extras.dat'.format(runs))
        self.rs = extras[3] * u.cm
        self.ms = extras[4] * u.g
        G = u.gravitational_constant
        rho = self.ms / (self.rs**3)
        self.tdyn = 1 / np.sqrt(G * rho)
        print('Tdyn = {:.2f}'.format(self.tdyn))

    def orbit(self, f=None, ax=None):
        if f == None:
            f, ax = plt.subplots(figsize=(12, 8))
        x = ((self.star[:, posx] - self.bh[:, posx]) * u.cm).in_units('AU')
        y = ((self.star[:, posy] - self.bh[:, posy]) * u.cm).in_units('AU')
        vel2 = 0
        for axis in [velx, vely, velz]:
            vel2 += ((self.star[:, axis] - self.bh[:, axis])
                     * u.cm / u.s)**2
        r2 = 0
        for axis in [posx, posy, posz]:
            r2 += ((self.star[:, axis] - self.bh[:, axis])
                   * u.cm)**2
        gpot = u.gravitational_constant * self.bh[:, mass] * u.gram
        E = -gpot / np.sqrt(r2) + 0.5 * vel2
        P = np.pi / np.sqrt(2) * gpot / (-E)**(3 / 2)
        print(np.mean(P).in_units('day'), np.std(P, ddof=1).in_units('day'))
        a = -gpot / 2 / E
        rp = np.sqrt(r2.min())
        print(rp)
        e = 1 - rp / a
        print(np.mean(e), np.std(e, ddof=1))
        pmean = np.mean(P).in_units('day').v
        plt.plot(x, y, label='{:.0f}'.format(pmean))

        '''ntdyn = 1
            arg = np.array([], dtype='i4')
            if self.tdyn > 0:
                while self.tdyn * ntdyn < self.T[-1]:
                    arg = np.append(arg, np.argwhere(
                        self.T > self.tdyn * ntdyn).flatten()[0])
                    ntdyn += 1
                plt.scatter(x[arg], y[arg], marker='+', s=50, color='k')'''
        ax.set_xlabel('X (AU)', fontsize=20)
        ax.set_ylabel('Y (AU)', fontsize=20)
        plt.axis('equal')
        f.tight_layout()

    def orbit_energy(self, f=None, ax=None):
        if f == None:
            f, ax = plt.figure(2, 1, figsize=(8, 16), sharex=True)
        x = ((self.star[:, posx] - self.bh[:, posx]) * u.cm).in_units('AU')
        y = ((self.star[:, posy] - self.bh[:, posy]) * u.cm).in_units('AU')
        vel2 = 0
        for axis in [velx, vely, velz]:
            vel2 += ((self.star[:, axis] - self.bh[:, axis])
                     * u.cm / u.s)**2
        r2 = 0
        for axis in [posx, posy, posz]:
            r2 += ((self.star[:, axis] - self.bh[:, axis])
                   * u.cm)**2
        gpot = u.gravitational_constant * self.bh[:, mass] * u.gram
        E = -gpot / np.sqrt(r2) + 0.5 * vel2
        P = np.pi / np.sqrt(2) * gpot / (-E)**(3 / 2)
        a = -gpot / 2 / E
        rp = np.sqrt(r2.min())
        e = 1 - rp / a
        pmean = np.mean(P).in_units('day').v
        ax[0].plot(np.sqrt(r2 / r2.min()), E,
                   label='{:.0f}'.format(pmean))
        ax[1].plot(np.sqrt(r2 / r2.min()),
                   (E - E[0]) / (gpot / np.sqrt(r2)).max())
        ax[1].set_xlabel('r (rp)', fontsize=20)
        ax[0].set_ylabel('E (erg/g)', fontsize=20)
        ax[1].set_ylabel('E (Ep)', fontsize=20)
        f.tight_layout()

    def delta_t(self):
        dT = self.T[1:] - self.T[:-1]
        # print(dT.min(), dT.max(), np.median(dT))
        f, ax = plt.subplots(2, 1, figsize=(12, 12))
        ax[0].hist(dT, bins=100)
        ax[1].plot(range(len(dT)), dT)
        if self.tdyn > 0:
            arg = np.argwhere(self.T > self.tdyn * 5).flatten()[0]
            ax[1].axvline(arg, color='k', linestyle='--')
        plt.tight_layout()
        plt.savefig('../../TDE_plot/{}.pdf'.format(self.runs))
        plt.show()


import argparse
parser = argparse.ArgumentParser(
    description='Download and analyze the evolution of sink particles.')
parser.add_argument('-d', dest='DIR',
                    help='Directory', default='ptgcliu')
parser.add_argument('--run', '-r', dest='runs',
                    help='Specific runs directory', default=['m0.8_p1_b0.6'], nargs='+')
parser.add_argument('--force', '-f', dest='force',
                    help='Force to rewrite pruned datafile', default=False, action='store_true')
parser.add_argument('--clean', '-c', dest='clean',
                    help='Clean pruned datafiles', default=False, action='store_true')
args = parser.parse_args()

#f1, ax1 = plt.subplots(figsize=(8, 8))
f2, ax2 = plt.subplots(2, 1, figsize=(8, 9))
for run in args.runs:
    test = SinkEvol(DIR=args.DIR, runs=run,
                    force=args.force, clean=args.clean)
    #test.orbit(f=f1, ax=ax1)
    test.orbit_energy(f=f2, ax=ax2)
# f1.legend()
plt.show()
