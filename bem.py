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

    def __init__(self, blade, u_inf, tsr, yaw):
        self.blade = blade
        self.u_inf = u_inf
        self.tsr = tsr
        self.yaw = np.radians(yaw)

    def get_prandtl(self, a, mu):
        '''
        Calculates the Prandtl loss factor.
        '''
        d = 2*np.pi/self.blade.n_blades*(1-a)/np.sqrt(self.tsr**2+1/(mu**2)*(1-a)**2)
        f_tip = 2/np.pi*np.arccos(np.exp(-np.pi*(1-mu)/d))
        mu_hub = self.blade.df['mu'][0]
        f_root = 2/np.pi*np.arccos(np.exp(-np.pi*(mu-mu_hub)/d))
        f = f_tip*f_root+1e-4
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
        fx = 0.5*w2_u2*self.u_inf**2*cx*self.blade.n_blades
        fx = fx*chord*dpsi/(2*np.pi)*dmu*self.blade.radius
        fy = 0.5*w2_u2*self.u_inf**2*cy*self.blade.n_blades
        fy = fy*chord*dpsi/(2*np.pi)*dmu*self.blade.radius
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
        return a, ap, fx, fy


class Rotor:
    '''
    Stores the BEM results over the rotor for the given conditions.
    '''

    def __init__(self, blade, n_r, n_az):
        '''
        Discretises the rotor.
        '''
        self.blade = blade
        self.n_r = n_r
        self.mu = np.linspace(0.2, 1, n_r)
        self.n_az = n_az
        self.psi = np.linspace(-180, 180, n_az)

    def solve(self, u_inf, tsr, yaw):
        '''
        Solves the BEM equations for the given conditions.
        '''
        solver = BEMSolver(self.blade, u_inf, tsr, yaw)
        a = np.zeros((self.n_r, self.n_az))
        ap = np.zeros((self.n_r, self.n_az))
        fx = np.zeros((self.n_r, self.n_az))
        fy = np.zeros((self.n_r, self.n_az))
        for i in range(self.n_r-1):
            for j in range(self.n_az-1):
                mu = 0.5*(self.mu[i]+self.mu[i+1])
                dmu = self.mu[i+1]-self.mu[i]
                psi = 0.5*(self.psi[j]+self.psi[j+1])
                dpsi = self.psi[j+1]-self.psi[j]
                solution = solver.solve(mu, dmu, psi, dpsi)
                a[i, j] = solution[0]
                ap[i, j] = solution[1]
                fx[i, j] = solution[2]
                fy[i, j] = solution[3]
            print(f'Position mu = {mu:.4f} has been solved.')
        self.a = a
        self.ap = ap
        self.fx = fx
        self.fy = fy


if __name__ == '__main__':
    blade = Blade()
    u_inf, tsr, yaw = 10, 6, 30
    # solver = BEMSolver(blade, u_inf, tsr, yaw)
    # a, ap, fx, fy = solver.solve(0.505, 0.005, 30, 20)
    rotor = Rotor(blade, 51, 21)
    rotor.solve(u_inf, tsr, yaw)
    # print(a, ap, fx, fy)
    # graph(blade, solver)
