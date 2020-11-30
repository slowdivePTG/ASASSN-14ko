import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt.units as u
import yt


class SinkEvol:
    def __init__(self, DIR, runs, force_download):
        self.user = DIR
        self.directory = runs

        import paramiko
        import os
        dirs = os.listdir('../../TDE_plot/')
        exist = False
        for i in dirs:
            temp = '{}_pruned_sinks_evol.dat'.format(runs)
            if temp in i:
                exist = True
                break
        if not exist or force_download:
            hostname = "fend01.hpc.ku.dk"
            port = 22
            username = "ptgcliu"
            password = "19980321lcPTG+"

            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(hostname, port, username, password, compress=True)
            sftp_client = client.open_sftp()

            file_cluster = '/storage/dark/{}/{}/sinks_evol.dat'.format(
                DIR, runs)
            filename = '../../TDE_plot/{}_sinks_evol.dat'.format(runs)
            sftp_client.get(file_cluster, filename)
            df = pd.read_table(filename, delim_whitespace=True)
            args = df[df['[01]time'] == '[01]time'].index
            cut = 0
            for i in range(len(args)):
                if eval(df['[01]time'][args[-i] + 1]) <= 1e-10:
                    cut = args[-i] + 1
                    break
            df = df[cut:]
            df = df[df['[01]time'] != '[01]time']
            df.to_csv('../../TDE_plot/{}_pruned_sinks_evol.dat'.format(runs),
                      sep=' ',
                      index=False,
                      header=True)

    def orbit(self):
        df = pd.read_table('../../TDE_plot/{}_pruned_sinks_evol.dat'.format(
            self.directory), delim_whitespace=True)
        T = np.array(df['[01]time'], dtype='f4')

        plt.figure(figsize=(12, 8))
        plt.scatter(
            df['[02]posx'][0::5000] * u.cm.in_units('AU'),
            df['[03]posy'][0::5000] * u.cm.in_units('AU'))
        plt.scatter(
            df['[02]posx'][1::500] * u.cm.in_units('AU'),
            df['[03]posy'][1::500] * u.cm.in_units('AU'),
            s=100**(np.log10(np.array(df['[01]time'][1::500], dtype='f4')) -
                    3.5))
        plt.xlabel('X (AU)', fontsize=20)
        plt.ylabel('Y (AU)', fontsize=20)
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig('{}.pdf'.format(self.directory))
        plt.show()


import argparse
parser = argparse.ArgumentParser(
    description='Download and analyze the evolution of sink particles.')
parser.add_argument('-d', dest='DIR',
                    help='Directory', default='ptgcliu')
parser.add_argument('--run', '-r', dest='runs',
                    help='Specific runs directory', default='M0.8_1')
parser.add_argument('--force', '-f', dest='force',
                    help='Specific runs directory', default=False, action='store_true')
args = parser.parse_args()

test = SinkEvol(DIR=args.DIR, runs=args.runs, force_download=args.force)
test.orbit()
