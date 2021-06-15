import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def data_binning(data, size=2, pri=False):  # size - day
    data_bin = []
    i = 0
    if pri:
        print(len(data))
    while i < len(data):
        j = i
        while j < len(data):
            if (data[j, 0] < data[i, 0] + size):
                j += 1
            else:
                break
        temp = data[i:j, :]
        if pri:
            print(i, j, temp[:, 0] - 2450000)
        if len(temp) > 1:
            arg = np.argwhere(
                abs(temp[:, 1] - temp[:, 1].mean()) <= 3 *
                np.std(temp[:, 1], ddof=1)).flatten()
            date_bin = temp[arg, 0].mean()
            weight_mag = (temp[arg, 1] / temp[arg, 2]**2) / \
                (temp[arg, 2]**(-2)).sum()
            weight_mag2 = (temp[arg, 1]**2 / temp[arg, 2]**2) / \
                (temp[arg, 2]**(-2)).sum()
            mag_bin = weight_mag.sum()
            magerr_bin0 = ((data[arg, 2]**(-2)).sum())**(-0.5)
            magerr_var = weight_mag2.sum() - mag_bin**2
            magerr_bin = (magerr_bin0**2 + magerr_var)**(0.5)
            i = j
        else:
            date_bin = data[i, 0]
            mag_bin = data[i, 1]
            magerr_bin = data[i, 2]
            i = i + 1
        data_bin.append([date_bin, mag_bin, magerr_bin])
    return np.array(data_bin)
