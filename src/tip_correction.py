import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import bem


plt.rcdefaults()
plt.rc('lines', markersize=4)
plt.rc('markers', fillstyle='none')
plt.rc('axes', grid=True)
plt.rc('legend', framealpha=0.7)
plt.rc('savefig', dpi=200, bbox='tight')
plt.close('all')


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
    return df


def tip_corr_results():
    '''
    Returns two rotor results (with and without Prandtl's tip loss correction).
    '''
    prandtl = read_csv('../results/tip-correction/prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
    no_prandtl = read_csv('../results/tip-correction/no-prandtl/u_inf-10-tsr-8-yaw-0-n_r-51-n_az-5_az_av.csv')
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


def print_results(prandtl, no_prandtl):
    '''
    Prints the results.
    '''
    print('Power coefficient CP')
    cp_wc = sum(prandtl['dcp'])
    cp_nc = sum(no_prandtl['dcp'])
    dcp = (cp_nc-cp_wc)/cp_wc
    print(f'Tip correction: {cp_wc:.4f}')
    print(f'No tip correction: {cp_nc:.4f}')
    print(f'Increase: {dcp*1e2:.2f} %')
    print('\nThrust coefficient CT')
    ct_wc = sum(prandtl['dct'])
    ct_nc = sum(no_prandtl['dct'])
    dct = (ct_nc-ct_wc)/ct_wc
    print(f'Tip correction: {ct_wc:.4f}')
    print(f'No tip correction: {ct_nc:.4f}')
    print(f'Increase: {dct*1e2:.2f} %')
    print('\nPower to thrust ratio CP/CT')
    cp_ct_wc = cp_wc/ct_wc
    cp_ct_nc = cp_nc/ct_nc
    dcp_ct = (cp_ct_nc-cp_ct_wc)/cp_ct_wc
    print(f'Tip correction: {cp_ct_wc:.4f}')
    print(f'No tip correction: {cp_ct_nc:.4f}')
    print(f'Increase: {dcp_ct*1e2:.2f} %')


if __name__ == '__main__':
    prandtl, no_prandtl = tip_corr_results()
    path = '../results/tip-correction'
    # plot_tip_corr(prandtl, no_prandtl, 'dcp_dmu', r'$dC_P/d\mu$ [-]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'dcs_dmu', r'$dC_S/d\mu$ [-]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'dct_dmu', r'$dC_T/d\mu$ [-]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'alpha', r'$\alpha$ [deg]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'phi', r'$\phi$ [deg]', path)
    plot_tip_corr(prandtl, no_prandtl, 'w2_u2', r'$w^2/u_{\infty}^2$ [deg]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'a', r'a [-]', path)
    # plot_tip_corr(prandtl, no_prandtl, 'f', r'Tip correction f [-]', path)
    print_results(prandtl, no_prandtl)
    plt.close('all')
