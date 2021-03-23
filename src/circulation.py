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


U_INF = 10
TSR = 8
R = 50
N_BLADES = 3


def read_csv(path):
    '''
    Returns the data frame of the given file.
    '''
    df = pd.read_csv(path)
    df['dct_dmu'] = df['dct']/df['dmu']
    df['ct_loc'] = df['dct_dmu']/(2*df['mu'])
    df['dcp_dmu'] = df['dcp']/df['dmu']
    df['dcs_dmu'] = df['dcp']/(TSR*df['mu']*df['dmu'])
    df['gamma_'] = df['gamma']*TSR/(U_INF*R)
    df['gamma_conv'] = N_BLADES*np.sqrt(df['w2_u2'])*df['gamma']*np.cos(df['phi'])/(U_INF*np.pi*df['mu']*R)
    return df


def plot_stagnation_enthalpy(df, field, label, path):
    '''
    Plots the stagnation entalpy over the blade.
    '''
    fig, ax = plt.subplots()
    ax.plot(df['mu'], df[field])
    ax.set_xlabel(r'$\mu$ [-]')
    ax.set_ylabel(label)
    ax.set_title(r'$\lambda = 8$ [-], $\gamma = 0$ [deg]')
    fig.show()
    if not os.path.isdir(path):
        os.makedirs(path)
    png_file = os.path.join(path, field+'.pdf')
    fig.savefig(png_file)


if __name__ == '__main__':
    plt.close('all')
    df = read_csv('../results/tip-correction/prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
    plot_stagnation_enthalpy(df, 'gamma_', r'$\lambda\Gamma/(u_{\infty} R)$ [-]', '../results/circulation')
    plot_stagnation_enthalpy(df, 'dct_dmu', r'$dC_T/d\mu$ [-]', '../results/circulation')
    plot_stagnation_enthalpy(df, 'gamma_conv', r'$w/u \Gamma /(uR)\cos \phi$ [-]', '../results/circulation')
    plot_stagnation_enthalpy(df, 'ct_loc', r'$C_T$ [-]', '../results/circulation')
