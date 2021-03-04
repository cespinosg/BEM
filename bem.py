import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
        d = 2*np.pi/self.blade.n_blades*(1-a)/np.sqrt(self.tsr**2+1/(mu**2)*(1-a)**2)
        f_tip = 2/np.pi*np.arccos(np.exp(-np.pi*(1-mu)/d))
        mu_hub = self.blade.df['mu'][0]
        f_root = 2/np.pi*np.arccos(np.exp(-np.pi*(mu-mu_hub)/d))
        f = f_tip*f_root+1e-4
        return f

    def get_flow_angle(self, a, ap, mu):
        '''
        Calculates the flow angle for the given induction factors.
        '''
        phi = np.arctan2(1-a, self.tsr*mu*(1+ap))
        return phi

    def solve(self, mu, tol=1e-9):
        '''
        Solves the BEM equations at the given section.
        '''
        chord, twist = self.blade.get_geo(mu)
        sigma_r = self.blade.n_blades*chord/(2*np.pi*mu*self.blade.radius)
        a0, ap0 = 1/3, 0
        fx, fy = 0, 0
        error = tol+1
        i = 0
        while tol < error:
            phi = self.get_flow_angle(a0, ap0, mu)
            alpha = np.degrees(phi)-twist-self.blade.pitch
            cl, cd = self.blade.get_cl_cd(alpha)
            w2_u2 = (1-a0)**2+(self.tsr*mu*(1+ap0))**2
            cx = cl*np.cos(phi)+cd*np.sin(phi)
            cy = cl*np.sin(phi)-cd*np.cos(phi)
            f = self.get_prandtl(a0, mu)
            a = w2_u2*sigma_r*cx/(4*(1-a0)*f)
            ap = w2_u2*sigma_r*cy/(4*(1-a0)*self.tsr*mu*f)
            error = max([abs(a-a0), abs(ap-ap0)])
            # Relaxing the induction factors is necessary for convergence near
            # the root
            a0, ap0 = 0.75*a0+0.25*a, 0.75*ap0+0.25*ap
            fx = 0.5*(w2_u2*100)*cx*chord #Uo2 = 100
            fy = 0.5*(w2_u2*100)*cy*chord #Uo2 = 100
            #print(i)
            i += 1
        return a, ap, fx, fy

def graph(blade, solver):
    Uo = 10 #m/s
    omega = Uo*solver.tsr/blade.radius
    dr = 0.01
    mu = np.arange(0.2+dr/2, 1.0-dr/2, dr) #midpoints of annuli
    print(len(mu))
    a, ap, fx, fy = np.zeros(len(mu)), np.zeros(len(mu)), np.zeros(len(mu)), np.zeros(len(mu))

    for i in range(len(mu)):
        a[i], ap[i], fx[i], fy[i] = solver.solve(mu[i])
        print(i)

    CT = np.sum(dr*fx*blade.n_blades/(0.5*Uo**2*np.pi*blade.radius**2))
    CP = np.sum(dr*fy*mu*blade.n_blades*blade.radius*omega/(0.5*Uo**3*np.pi*blade.radius**2))
    print("CT = ", CT)
    print("CP = ", CP)

    '''
    plot a and a'
    '''
    fig1 = plt.figure(figsize=(12, 6))
    plt.title('Axial and tangential induction')
    plt.plot(mu, a, 'b-', label=r'$a$')
    plt.plot(mu, ap, 'g--', label=r'$a^,$')
    plt.grid()
    plt.xlabel('r/R')
    plt.legend()
    plt.savefig("a_aip.png")
    plt.show()

    '''
    plot Fax and Faz
    '''
    fig1 = plt.figure(figsize=(12, 6))
    plt.title(r'Normal and tagential force, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$')
    plt.plot(mu, fx/(0.5*Uo**2*blade.radius), 'b-', label=r'Faxial')
    plt.plot(mu, fy/(0.5*Uo**2*blade.radius), 'g--', label=r'Fazimuthal')
    plt.grid()
    plt.xlabel('μ')
    plt.legend()
    plt.savefig("Fax_Faz.png") 
    plt.show()

if __name__ == '__main__':
    blade = Blade()
    solver = BEMSolver(blade, 6)
    #a, ap, fx, fy = solver.solve(0.505)
    #print(a, ap, fx, fy)
    graph(blade, solver)
