import numpy as np
import pandas as pd


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

    def __init__(self, blade, tsr):
        self.blade = blade
        self.tsr = tsr

    def get_prandtl(self, a, mu):
        '''
        Calculates the Prandtl loss factor.
        '''
        d = 2*np.pi/self.blade.n_blades*(1-a)/np.sqrt(self.tsr**2+(1-a)**2)
        f_tip = 2/np.pi*np.arccos(np.exp(-np.pi*(1-mu)/d))
        mu_hub = self.blade.df['mu'][0]
        f_root = 2/np.pi*np.arccos(np.exp(-np.pi*(mu-mu_hub)/d))
        f = f_tip*f_root
        return f

    def get_flow_angle(self, a, ap, mu):
        '''
        Calculates the flow angle for the given induction factors.
        '''
        f = self.get_prandtl(a, mu)
        phi = np.arctan2(1-a/f, self.tsr*mu*(1+ap/f))
        return phi

    def solve(self, mu, tol=1e-6):
        '''
        Solves the BEM equations at the given section.
        '''
        chord, twist = self.blade.get_geo(mu)
        sigma_r = self.blade.n_blades*chord/(2*np.pi*mu*self.blade.radius)
        a0, ap0 = 1/3, 0
        error = tol+1
        i = 0
        while tol < error:
            phi = self.get_flow_angle(a0, ap0, mu)
            alpha = np.degrees(phi)-twist-self.blade.pitch
            cl, cd = self.blade.get_cl_cd(alpha)
            w2_u2 = (1-a0)**2+(self.tsr*mu*(1+ap0))**2
            cx = cl*np.cos(phi)+cd*np.sin(phi)
            cy = cl*np.sin(phi)-cd*np.cos(phi)
            a = w2_u2*sigma_r*cx/(4*(1-a0))
            ap = w2_u2*sigma_r*cy/(4*(1-a0)*self.tsr*mu)
            error = max([abs(a-a0), abs(ap-ap0)])
            a0, ap0 = a, ap
            # print(i)
            # i += 1
        return a, ap

if __name__ == '__main__':
    blade = Blade()
    solver = BEMSolver(blade, 10)
    a, ap = solver.solve(0.25)
