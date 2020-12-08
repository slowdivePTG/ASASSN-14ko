import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt.units as u
import yt
import os


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
            temp = '{}_extras.dat'.format(runs)
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
                      header=True)
        df = pd.read_table(self.output, delim_whitespace=True)

        self.star = df[1::2]
        self.bh = df[0::2]
        self.T = np.array(self.star['[01]time'], dtype='f4')
        os.system('rm ../../TDE_plot/sinks_evol.dat')
        self.runs = runs

        extras = np.loadtxt('../../TDE_plot/{}_extras.dat'.format(runs))
        self.rs = extras[3] * u.cm
        self.ms = extras[4] * u.g
        G = u.gravitational_constant
        rho = self.ms / (self.rs**3)
        self.tdyn = 1 / np.sqrt(G * rho)
        print('Tdyn = {:.2f} s'.format(self.tdyn))

    def orbit(self):
        plt.figure(figsize=(12, 8))
        plt.scatter(
            self.bh['[02]posx'][::200] * u.cm.in_units('AU'),
            self.bh['[03]posy'][::200] * u.cm.in_units('AU'))
        plt.scatter(
            self.star['[02]posx'][::200] * u.cm.in_units('AU'),
            self.star['[03]posy'][::200] * u.cm.in_units('AU'),
            s=self.star['[01]time'][::200]**2 * 1e-7)
        ntdyn = 1
        arg = np.array([], dtype='i4')
        if self.tdyn > 0:
            while self.tdyn * ntdyn < self.T[-1]:
                arg = np.append(arg, np.argwhere(
                    self.T > self.tdyn * ntdyn).flatten()[0])
                ntdyn += 1
            relaxx = np.array(self.star['[02]posx'])[arg] * u.cm.in_units('AU')
            relaxy = np.array(self.star['[03]posy'])[arg] * u.cm.in_units('AU')
            plt.scatter(relaxx, relaxy, marker='+', s=50, color='k')
        plt.xlabel('X (AU)', fontsize=20)
        plt.ylabel('Y (AU)', fontsize=20)
        plt.axis('equal')
        plt.tight_layout()
        plt.show()

    def delta_t(self):
        dT = self.T[1:] - self.T[:-1]
        #print(dT.min(), dT.max(), np.median(dT))
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
                    help='Specific runs directory', default='m0.8_p1_b0.6')
parser.add_argument('--force', '-f', dest='force',
                    help='Force to rewrite pruned datafile', default=False, action='store_true')
parser.add_argument('--clean', '-c', dest='clean',
                    help='Clean pruned datafiles', default=False, action='store_true')
args = parser.parse_args()

test = SinkEvol(DIR=args.DIR, runs=args.runs,
                force=args.force, clean=args.clean)
test.orbit()
test.delta_t()
