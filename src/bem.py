import os

import numpy as np
import pandas as pd
# from graph import *


class Blade():
    '''
    Stores the blade geometry distribtuion over the span and airfoil polar.
    '''

    def __init__(self):
        self.csv = 'blade-geo.csv'
        self.radius = 50
        self.n_blades = 3
        self.pitch = -2
        df = pd.DataFrame()
        df['mu'] = np.linspace(0.2, 1)
        df['z'] = self.radius*df['mu']
        df['twist'] = 14*(1-df['mu'])
        df['chord'] = 3*(1-df['mu'])+1
        self.df = df
        df = df.set_index('mu')
        df.to_csv(self.csv)
        self.polar = pd.read_csv('polar-DU95W180.csv')

    def get_geo(self, mu):
        '''
        Returns the blade geometry at the given mu.
        '''
        chord = np.interp(mu, self.df['mu'], self.df['chord'])
        twist = np.interp(mu, self.df['mu'], self.df['twist'])
        return chord, twist

    def get_cl_cd(self, alpha):
        '''
        Returns the aerodynamic coefficients for the given AoA.
        '''
        cl = np.interp(alpha, self.polar['alpha'], self.polar['cl'])
        cd = np.interp(alpha, self.polar['alpha'], self.polar['cd'])
        return cl, cd


class BEMSolver():
    '''
    Solves the BEM equations for the given geometry and conditions.
    '''

    keys = ['a', 'ap', 'f', 'cx', 'cy', 'fx', 'fy', 'phi', 'w2_u2', 'alpha',
            'cl', 'cd', 'gamma', 'dmu', 'dpsi']

    def __init__(self, blade, u_inf, tsr, yaw, prandtl=True):
        self.blade = blade
        self.u_inf = u_inf
        self.tsr = tsr
        self.yaw = np.radians(yaw)
        self.prandtl = prandtl
        self.name = f'u_inf-{u_inf}-tsr-{tsr}-yaw-{yaw}'

    def get_prandtl(self, a, mu):
        '''
        Calculates the Prandtl loss factor.
        '''
        if self.prandtl:
            d = 2*np.pi/self.blade.n_blades*(1-a)/np.sqrt(self.tsr**2+1/(mu**2)*(1-a)**2)
            f_tip = 2/np.pi*np.arccos(np.exp(-np.pi*(1-mu)/d))
            mu_hub = self.blade.df['mu'][0]
            f_root = 2/np.pi*np.arccos(np.exp(-np.pi*(mu-mu_hub)/d))
            f = f_tip*f_root+1e-4
        else:
            f = 1
        return f

    def get_flow_angle(self, u_a, u_t):
        '''
        Calculates the flow angle for the given the axial and tangential
        velocities.
        '''
        phi = np.arctan2(u_a, u_t)
        return phi

    def get_velocities(self, a, ap, mu, psi):
        '''
        Returns the axial and tangential velocities.
        '''
        chi = (0.6*a+1)*self.yaw
        k = 2*np.tan(0.5*chi)
        u_a = np.cos(self.yaw)-a*(1+k*mu*np.sin(psi))
        u_t = self.tsr*mu*(1+ap)-np.sin(self.yaw)*np.cos(psi)
        return u_a, u_t

    def get_forces(self, chord, twist, a0, ap0, mu, psi, dmu, dpsi):
        '''
        Calculates the forces at the given position.
        '''
        u_a, u_t = self.get_velocities(a0, ap0, mu, psi)
        phi = self.get_flow_angle(u_a, u_t)
        alpha = np.degrees(phi)-twist-self.blade.pitch
        cl, cd = self.blade.get_cl_cd(alpha)
        w2_u2 = u_a**2+u_t**2
        cx = cl*np.cos(phi)+cd*np.sin(phi)
        cy = cl*np.sin(phi)-cd*np.cos(phi)
        fx = 0.5*w2_u2*self.u_inf**2*cx
        fx = fx*self.blade.n_blades*chord*dpsi/(2*np.pi)*dmu*self.blade.radius
        fy = 0.5*w2_u2*self.u_inf**2*cy
        fy = fy*self.blade.n_blades*chord*dpsi/(2*np.pi)*dmu*self.blade.radius
        return fx, fy

    def get_a(self, CT):
        '''
        Calculates a given CT applies glauert correction.
        '''
        CT1 = 1.816
        CT2 = 2*np.sqrt(CT1) - CT1
        if(CT < CT2):
            a = 0.5 - 0.5*np.sqrt(1-CT)
        elif(CT >= CT2):
            a = 1 + (CT-CT1)/(4*np.sqrt(CT1)-4)
        return a

    def solve(self, mu, dmu, psi, dpsi, tol=1e-9):
        '''
        Solves the BEM equations at the given section.
        '''
        psi = np.radians(psi)
        dpsi = np.radians(dpsi)
        chord, twist = self.blade.get_geo(mu)
        a0, ap0 = 1/3, 0
        area = np.pi*(((mu+0.5*dmu)*self.blade.radius)**2
                      - ((mu-0.5*dmu)*self.blade.radius)**2)
        area = area*dpsi/(2*np.pi)
        error = tol+1
        i = 0
        while tol < error:
            f = self.get_prandtl(a0, mu)
            fx, fy = self.get_forces(chord, twist, a0, ap0, mu, psi, dmu, dpsi)
            CT = fx/(0.5*area*self.u_inf**2)
            cfy = fy/(0.5*area*self.u_inf**2)
            a = self.get_a(CT)
            ap = cfy/(4*(1-a)*self.tsr*mu)
            a = a/f
            ap = ap/f
            error = max([abs(a-a0), abs(ap-ap0)])
            # Relaxing the induction factors is necessary for convergence near
            # the root
            a0, ap0 = 0.75*a0+0.25*a, 0.75*ap0+0.25*ap
            # print(i)
            i += 1
        f = self.get_prandtl(a, mu)
        u_a, u_t = self.get_velocities(a, ap, mu, psi)
        phi = self.get_flow_angle(u_a, u_t)
        alpha = np.degrees(phi)-twist-self.blade.pitch
        cl, cd = self.blade.get_cl_cd(alpha)
        cx = cl*np.cos(phi)+cd*np.sin(phi)
        cy = cl*np.sin(phi)-cd*np.cos(phi)
        w2_u2 = u_a**2+u_t**2
        gamma = 0.5*self.u_inf*chord*cl
        solution = [a, ap, f, cx, cy, fx, fy, phi, w2_u2, alpha, cl, cd,
                    gamma, dmu, dpsi]
        solution = dict(zip(self.keys, solution))
        return solution


class Rotor:
    '''
    Stores the BEM results over the rotor for the given conditions.
    '''

    def __init__(self, n_r, n_az):
        '''
        Discretises the rotor with the given number of nodes.
        '''
        self.name = f'n_r-{n_r}-n_az-{n_az}'
        self.n_r = n_r
        self.mu_n = np.linspace(0.2, 1, n_r)  # nodes
        self.mu_e = 0.5*(self.mu_n[1:]+self.mu_n[:-1])  # elements
        self.dmu = self.mu_n[1:]-self.mu_n[:-1]
        self.n_az = n_az
        self.psi_n = np.linspace(-180, 180, n_az)  # nodes
        self.psi_e = 0.5*(self.psi_n[1:]+self.psi_n[:-1])  # elements
        self.dpsi = self.psi_n[1:]-self.psi_n[:-1]

    def solve(self, solver):
        '''
        Solves the BEM equations with the given solver.
        '''
        self.solver = solver
        index = pd.MultiIndex.from_product([self.mu_e, self.psi_e],
                                           names=['mu', 'azimuth'])
        df = pd.DataFrame(columns=solver.keys, index=index)
        for i in range(self.n_r-1):
            mu = self.mu_e[i]
            dmu = self.dmu[i]
            for j in range(self.n_az-1):
                psi = self.psi_e[j]
                dpsi = self.dpsi[j]
                df.loc[(mu, psi)] = solver.solve(mu, dmu, psi, dpsi)
            print(f'Position mu = {mu:.4f} has been solved.')
        self.df = df
        self.az_av = self.azimuth_average()

    def azimuth_average(self):
        '''
        Returns the fields averaged over the azimuth.
        '''
        cols = self.df.columns
        df = pd.DataFrame(columns=cols, index=self.mu_e)
        df.index.name = 'mu'
        for mu in self.mu_e:
            av = {c: sum(self.df.loc[mu][c]*self.dpsi)/sum(self.dpsi) for c in cols}
            df.loc[mu] = av
        df['fx'] = df['fx']*(self.n_az-1)
        df['fy'] = df['fy']*(self.n_az-1)
        return df

    def to_csv(self, path):
        '''
        Saves the results in the given directory.
        '''
        if not os.path.isdir(path):
            os.makedirs(path)
        csv_file = self.solver.name+'-'+self.name
        self.df.to_csv(os.path.join(path, csv_file+'.csv'))
        self.az_av.to_csv(os.path.join(path, csv_file+'_az_av.csv'))


if __name__ == '__main__':
    blade = Blade()
    u_inf, tsr, yaw = 10, 8, 0
    solver = BEMSolver(blade, u_inf, tsr, yaw, prandtl=False)
    # solution = solver.solve(0.505, 0.005, 30, 20)
    # print(solution)
    rotor = Rotor(51, 5)
    rotor.solve(solver)
    # check the azimuth average
    # av = rotor.azimuth_average()
    # az0 = list(set(rotor.df.index.get_level_values('azimuth')))[0]
    # print(max(abs(av['a']-rotor.df['a'].xs(az0, level='azimuth'))))
    rotor.to_csv('../results/tip-correction/no-prandtl')
    # graph(blade, solver)
