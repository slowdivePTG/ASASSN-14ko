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
            temp = '{}_flash.par'.format(runs)
            if temp in i:
                exist = True
                break

        if not exist or force:

            hostname = "fend02.hpc.ku.dk"
            port = 22
            username = "ptgcliu"
            password = "19980321lcPTG+"

            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(hostname, port, username, password, compress=True)
            sftp_client = client.open_sftp()
            sftp_client.get(
                '/storage/dark/ptgcliu/{}/flash.par'.format(self.runs), '../../TDE_plot/{}_flash.par'.format(runs))

            file_cluster = '/storage/dark/{}/{}/sinks_evol.dat'.format(
                DIR, self.runs)
            # print(file_cluster)
            sftp_client.get(file_cluster, filename)
            df = pd.read_table(filename, delim_whitespace=True)
            args = df[df['[01]time'] == '[01]time'].index
            cut = 0
            while len(args) != 0:
                sub = df['[01]time'][args[-1] + 1]
                df_before = df[df.index < args[-1]]
                arg_temp = df_before[df_before['[01]time'] == sub].index
                if len(arg_temp) > 0:
                    remove = range(arg_temp[0], args[-1] + 1)
                else:
                    remove = [args[-1]]
                df = df.drop(remove)
                args = df[df['[01]time'] == '[01]time'].index
            df = df[cut:]
            df = df[df['[01]time'] != '[01]time']
            df.to_csv(self.output,
                      sep=' ',
                      index=False,
                      header=True)
        else:
            df = pd.read_table(self.output, delim_whitespace=True)

        self.star = df[1::2]
        self.bh = df[0::2]
        self.T = np.array(self.star['[01]time'], dtype='f4')
        os.system('rm ../../TDE_plot/sinks_evol.dat')
        self.runs = runs

        with open('../../TDE_plot/{}_flash.par'.format(runs)) as f:
            while True:
                l = f.readline()
                if 'tdyn' in l:
                    cut = l.find('tdyn = ') + len('tdyn = ')
                    try:
                        self.tdyn = eval(l[cut:])
                    except:
                        s = l.find(' s')
                        self.tdyn = eval(l[cut:s])
                    print('Tdyn = {} s'.format(self.tdyn))
                    break

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
        print(dT.min(), dT.max(), np.median(dT))
        f, ax = plt.subplots(2, 1, figsize=(12, 12))
        ax[0].hist(dT, bins=100)
        ax[1].plot(range(len(dT)), dT)
        ax[1].axvline(self.tdyn * 5, color='k', linestyle='--')
        plt.tight_layout()
        # plt.savefig('../../TDE_plot{}.pdf'.format(self.runs))
        plt.show()


import argparse
parser = argparse.ArgumentParser(
    description='Download and analyze the evolution of sink particles.')
parser.add_argument('-d', dest='DIR',
                    help='Directory', default='ptgcliu')
parser.add_argument('--run', '-r', dest='runs',
                    help='Specific runs directory', default='M0.8_1')
parser.add_argument('--force', '-f', dest='force',
                    help='Force to rewrite pruned datafile', default=False, action='store_true')
parser.add_argument('--clean', '-c', dest='clean',
                    help='Clean pruned datafiles', default=False, action='store_true')
args = parser.parse_args()

test = SinkEvol(DIR=args.DIR, runs=args.runs,
                force=args.force, clean=args.clean)
#test.orbit()
test.delta_t()
