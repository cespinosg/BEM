import os

import matplotlib.pyplot as plt
import pandas as pd

import bem


plt.rcdefaults()
plt.rc('lines', markersize=4)
plt.rc('markers', fillstyle='none')
plt.rc('axes', grid=True)
plt.rc('legend', framealpha=0.7)
plt.rc('savefig', dpi=200, bbox='tight')
plt.close('all')


def tip_corr_results():
    '''
    Returns two rotor results (with and without Prandtl's tip loss correction).
    '''
    prandtl = pd.read_csv('../results/tip-correction/prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
    no_prandtl = pd.read_csv('../results/tip-correction/no-prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
    return prandtl, no_prandtl


def plot_tip_corr(prandtl, no_prandtl, field, label, path):
    '''
    Compares the results with and without tip correction.
    '''
    fig, ax = plt.subplots()
    ax.plot(prandtl['mu'], prandtl[field], label='Tip correction')
    ax.plot(no_prandtl['mu'], no_prandtl[field],
               label='No tip correction')
    ax.set_ylabel(label)
    ax.set_xlabel(r'$\mu$ [-]')
    ax.legend()
    ax.set_title(r'$\lambda=8$ [-], $\gamma=0$ [deg]')
    fig.show()
    if not os.path.isdir(path):
        os.makedirs(path)
    png_file = os.path.join(path, field+'.pdf')
    fig.savefig(png_file)


if __name__ == '__main__':
    prandtl, no_prandtl = tip_corr_results()
    path = '../results/tip-correction'
    plot_tip_corr(prandtl, no_prandtl, 'fy', r'$f_y$ [N/(kg/m3)]', path)
    plot_tip_corr(prandtl, no_prandtl, 'alpha', r'$\alpha$ [deg]', path)
    plot_tip_corr(prandtl, no_prandtl, 'a', r'a [-]', path)
    plot_tip_corr(prandtl, no_prandtl, 'f', r'Tip correction f [-]', path)
    plt.close('all')
