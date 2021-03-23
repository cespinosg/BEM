import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


plt.rcdefaults()
plt.rc('lines', markersize=4)
plt.rc('markers', fillstyle='none')
plt.rc('axes', grid=True)
plt.rc('legend', framealpha=0.7)
plt.rc('savefig', dpi=200, bbox='tight')


TSR = 8


def read_csv(path):
    '''
    Returns the data frame of the given file.
    '''
    df = pd.read_csv(path)
    df['dct_dmu'] = df['dct']/df['dmu']
    df['dcp_dmu'] = df['dcp']/df['dmu']
    df['dcs_dmu'] = df['dcp']/(TSR*df['mu']*df['dmu'])
    df['phi'] = np.degrees(df['phi'])
    df['pt_u'] = 0*df['mu']+1
    df['pt_d'] = (1-2*df['a'])**2
    return df


def plot_stagnation_enthalpy(df, path):
    '''
    Plots the stagnation entalpy over the blade.
    '''
    fig, ax = plt.subplots()
    ax.plot(df['mu'], df['pt_u'], label='Up-stream')
    ax.plot(df['mu'], df['pt_d'], label='Down-stream')
    ax.set_xlabel(r'$\mu$ [-]')
    ax.set_ylabel(r'$2(p_t-p_a)/u_{\infty}^2$ [-]')
    ax.set_title(r'$\lambda = 8$ [-], $\gamma = 0$ [deg]')
    ax.legend()
    fig.show()
    if not os.path.isdir(path):
        os.makedirs(path)
    png_file = os.path.join(path, 'stagnation-enthalpy.pdf')
    fig.savefig(png_file)


if __name__ == '__main__':
    plt.close('all')
    df = read_csv('../results/tip-correction/prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
    plot_stagnation_enthalpy(df, '../results')
