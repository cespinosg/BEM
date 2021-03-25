import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import sys

plt.rcdefaults()
plt.rc('lines', markersize=4)
plt.rc('markers', fillstyle='none')
plt.rc('axes', grid=True)
plt.rc('legend', framealpha=0.7)
plt.rc('savefig', dpi=200, bbox='tight')
plt.close('all')

def graph(data6, data8, data10, tsr, yaw, n_az, path):#blade, solver):
    Uo = 10 #m/s
    radius = 50
    n_blades = 3
    omega6 = Uo*6/radius
    omega8 = Uo*8/radius
    omega10 = Uo*10/radius
    #dmu = 0.01
    #mu = np.arange(0.2+dmu/2, 1.0, dmu) #midpoints of annuli
    #print(len(mu))
    #a, ap, fx, fy = np.zeros(len(mu)), np.zeros(len(mu)), np.zeros(len(mu)), np.zeros(len(mu))

    #for i in range(len(mu)):
    #    a[i], ap[i], fx[i], fy[i] = solver.solve(mu[i], dmu)
    #    print(i)


    airfoil = 'polar-DU95W180.csv'
    data1=pd.read_csv(airfoil, header=0,
                        names = ["alpha", "cl", "cd", "cm"])
    polar_alpha = data1['alpha'][:]
    polar_cl = data1['cl'][:]
    polar_cd = data1['cd'][:]
    fig1 = plt.figure()
    plt.plot(polar_alpha, polar_cl, label=r'$C_l$')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$C_l$')
    plt.legend()
    plt.title(rf'$C_l\ vs\ \alpha$')
    plt.savefig(path+'Cl_v_alpha'+'.pdf')
    #plt.show()
    #plt.close('all')
    fig1 = plt.figure()
    plt.plot(polar_alpha, polar_cd, label=r'$C_d$')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$C_d$')
    plt.legend()
    plt.title(rf'$C_d\ vs\ \alpha$')
    plt.savefig(path+'Cd_v_alpha'+'.pdf')
    #plt.show()
    fig1 = plt.figure()
    plt.plot(polar_cd, polar_cl, label=r'$C_d$')
    plt.xlabel(r'$C_d$')
    plt.ylabel(r'$C_l$')
    plt.legend()
    plt.title(rf'$C_l\ vs\ C_d$')
    plt.savefig(path+'Cd_v_Cl'+'.pdf')
    #plt.show()


    mu6 = data6['mu']
    dmu6 = data6['dmu']
    a6 = data6['a']
    ap6 = data6['ap']
    alpha6 = data6['alpha']
    phi6 = data6['phi']
    cx6 = data6['cx']
    cy6 = data6['cy']
    fx6 = data6['fx']
    fy6 = data6['fy']

    mu8 = data8['mu']
    dmu8 = data8['dmu']
    a8 = data8['a']
    ap8 = data8['ap']
    alpha8 = data8['alpha']
    phi8 = data8['phi']
    cx8 = data8['cx']
    cy8 = data8['cy']
    fx8 = data8['fx']
    fy8 = data8['fy']

    mu10 = data10['mu']
    dmu10 = data10['dmu']
    a10 = data10['a']
    ap10 = data10['ap']
    alpha10 = data10['alpha']
    phi10 = data10['phi']
    cx10 = data10['cx']
    cy10 = data10['cy']
    fx10 = data10['fx']
    fy10 = data10['fy']

    azimuth = 0;
    if( yaw!=0 & n_az == 51 ):
        azimuth = data['azimuth']
        azimuth = azimuth[:n_az-1]

    dr6 = dmu6*radius#((mu+dmu/2) - (mu-dmu/2))*blade.radius
    CT6 = np.sum(fx6/(0.5*(Uo**2)*np.pi*radius**2))
    CP6 = np.sum(fy6*mu6*radius*omega6/(0.5*(Uo**3)*np.pi*radius**2))
    print("CT6 = ", CT6)
    print("CP6 = ", CP6)

    dr8 = dmu8*radius#((mu+dmu/2) - (mu-dmu/2))*blade.radius
    CT8 = np.sum(fx8/(0.5*(Uo**2)*np.pi*radius**2))
    CP8 = np.sum(fy8*mu8*radius*omega8/(0.5*(Uo**3)*np.pi*radius**2))
    print("CT8 = ", CT8)
    print("CP8 = ", CP8)

    dr10 = dmu10*radius#((mu+dmu/2) - (mu-dmu/2))*blade.radius
    CT10 = np.sum(fx10/(0.5*(Uo**2)*np.pi*radius**2))
    CP10 = np.sum(fy10*mu10*radius*omega10/(0.5*(Uo**3)*np.pi*radius**2))
    print("CT10 = ", CT10)
    print("CP10 = ", CP10)

    '''
    plot alpha and phi
    '''
    fig1 = plt.figure()
    #plt.title('Angle of attack and Inflow angle')
    plt.plot(mu6, alpha6, label=r'$angle\ of\ attack, tsr=6$')
    plt.plot(mu8, alpha8, label=r'$angle\ of\ attack, tsr=8$')
    plt.plot(mu10, alpha10, label=r'$angle\ of\ attack, tsr=10$')
    plt.plot(mu6, phi6*180/np.pi, label=r'$inflow\ angle, tsr=6$')
    plt.plot(mu8, phi8*180/np.pi, label=r'$inflow\ angle, tsr=8$')
    plt.plot(mu10, phi10*180/np.pi, label=r'$inflow\ angle, tsr=10$')
    #plt.grid()
    plt.xlabel(r'$\mu$ [-]')
    plt.ylabel('angle [deg]')
    plt.legend()
    plt.title(rf'$\lambda={6,8,10}$ [-], $\gamma={yaw}$ [deg]')
    plt.savefig(path+'alpha_phi'+f'-yaw_{yaw}'+'.pdf')
    #plt.show()


    '''
    plot a and a'
    '''
    fig1 = plt.figure()#figsize=(12, 6))
    plt.title('Axial and tangential induction')
    plt.plot(mu6, a6, label=r'$a, tsr=6$')
    plt.plot(mu8, a8, label=r'$a, tsr=8$')
    plt.plot(mu10, a10, label=r'$a, tsr=10$')
    plt.plot(mu6, ap6, label=r'$a^,, tsr=6$')
    plt.plot(mu8, ap8, label=r'$a^,, tsr=8$')
    plt.plot(mu10, ap10, label=r'$a^,, tsr=10$')
    #plt.grid()
    plt.xlabel(r'$\mu$ [-]')
    plt.legend()
    plt.title(rf'$\lambda={6,8,10}$ [-], $\gamma={yaw}$ [deg]')
    plt.savefig(path+'a_aip'+f'-yaw_{yaw}'+'.pdf')
    #plt.show()

    '''
    plot Fax and Faz
    '''
    fig1 = plt.figure()#figsize=(12, 6))
    plt.title(r'Normal and tagential force, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$')
    plt.plot(mu6, fx6/(dr6*3*0.5*Uo**2*radius),  label=r'Faxial, tsr=6')
    plt.plot(mu8, fx8/(dr8*3*0.5*Uo**2*radius),  label=r'Faxial, tsr=8')
    plt.plot(mu10, fx10/(dr10*3*0.5*Uo**2*radius),  label=r'Faxial, tsr=10')
    plt.plot(mu6, fy6/(dr6*3*0.5*Uo**2*radius),  label=r'Fazimuthal, tsr=6')
    plt.plot(mu8, fy8/(dr8*3*0.5*Uo**2*radius),  label=r'Fazimuthal, tsr=8')
    plt.plot(mu10, fy10/(dr10*3*0.5*Uo**2*radius),  label=r'Fazimuthal, tsr=10')
    #plt.grid()
    plt.xlabel(r'$\mu$ [-]')
    plt.legend()
    plt.title(rf'$\lambda={6,8,10}$ [-], $\gamma={yaw}$ [deg]')
    plt.savefig(path+'Fax_Faz'+f'-yaw_{yaw}'+'.pdf')
    #plt.show()

    plt.close('all')

    if(yaw!=0 & n_az == 51):
        '''
        plot a on contour plot
        '''
        #r = mu, ang = azimuth, value = a
        MU, AZIMUTH = np.meshgrid(mu[::n_az-1], azimuth[:n_az-1])
        X, Y = MU*np.sin(AZIMUTH*np.pi/180), MU*np.cos(AZIMUTH*np.pi/180)
        x, y = mu[::n_az-1]*np.sin(azimuth[:n_az-1]*np.pi/180), mu[::n_az-1]*np.cos(azimuth[:n_az-1]*np.pi/180)

        a = np.array(a)
        anew = a[::n_az-1].copy()
        for i in range(1,n_az-1):
            #print(i)
            anew = np.concatenate([anew, a[i::n_az-1].copy()])
        anew = np.reshape( anew, (n_az-1, int(len(anew)/(n_az-1))) )
        fig1 = plt.figure()
        con = plt.contourf(X,Y,anew)
        col = plt.colorbar()
        col.ax.set_ylabel('   a [-]',rotation=0)

        plt.title(rf'$\lambda={tsr}$ [-], $\gamma={yaw}$ [deg]')
        #plt.savefig(path+'contour.png')
        plt.savefig(path+'contour-'+f'n_az-{n_az}'+'.pdf')
        plt.show()



if __name__ == '__main__':
    #print(sys.argv)
    uinf = 10
    tsr = 6
    yaw = 30
    n_r = 51
    n_az = 5
    file6 = f'../results/yaw/u_inf-{uinf}-tsr-{tsr}-yaw-{yaw}-n_r-{n_r}-n_az-{n_az}_az_av.csv'
    tsr = 8
    file8 = f'../results/yaw/u_inf-{uinf}-tsr-{tsr}-yaw-{yaw}-n_r-{n_r}-n_az-{n_az}_az_av.csv'
    tsr = 10
    file10 = f'../results/yaw/u_inf-{uinf}-tsr-{tsr}-yaw-{yaw}-n_r-{n_r}-n_az-{n_az}_az_av.csv'
    if(yaw!=0):
        file = f'../results/yaw/u_inf-{uinf}-tsr-{tsr}-yaw-{yaw}-n_r-{n_r}-n_az-{n_az}.csv'
    path = '../results/yaw/'
    #file = sys.argv[1]
    #print(file)
    data6 = pd.read_csv(file6)
    data8 = pd.read_csv(file8)
    data10 = pd.read_csv(file10)
    graph(data6, data8, data10, tsr, yaw, n_az, path)
