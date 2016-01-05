import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.constants
import Gnuplot as gp


mm = 0.001
g = q = 2. * mm 
p = 4.*mm
eps0 = scipy.constants.epsilon_0
eps1 = eps3 = 10.*eps0
eps2 = eps0
e0 = scipy.constants.e

def SpaceChargeField(r,phi,z,rp,phip,zp):
    #eps0 = 8.854e-12
    #eps1 = 10.*eps0 #bakelite
    #eps2 = eps0
    #eps3 = eps1
    #e0 = 1.60217657e-19
    P2 = r*r - 2*r*rp*np.cos(phi-phip) + rp*rp
    #g = 0.2 #cm
    
    E = ( e0/(4*np.pi*eps2) ) * (   ( (z-zp)/((P2 + (z-zp)*(z-zp))**1.5) )  -  ( ((eps2-eps3)/(eps2+eps3))*(2*g-z-zp)/((P2+(2*g-z-zp)**2)**1.5) )     )
    
    return E

def Plot_field():
    r = pl.linspace(-1.*mm,1.*mm,50)
    z = pl.linspace(0,g,50)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.gca().patch.set_facecolor('white')
    ax.w_xaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_yaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_zaxis.set_pane_color((1., 1., 1., 1.0))
    #x = y = np.arange(-3.0, 3.0, 0.05)
    X, Y = np.meshgrid(z, r)
    zs = np.array([SpaceChargeField(r,0,z,0,0,0.5*mm) for z,r in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)
    
    ax.set_xlabel(r'$z\,[mm]$', fontsize=20)
    ax.set_ylabel(r'$r\,[mm]$', fontsize=20)
    ax.set_zlabel(r'$E\,[V/m]$', fontsize=20)
    
    ax.plot_wireframe(X*1000, Y*1000, Z*0.001,rstride=1,cstride=1)
    #plt.savefig("field_simplified.pdf")
    plt.show()
    
    
    
def Plot_field_gp():
    r = pl.linspace(-1.*mm,1.*mm,50)
    z = pl.linspace(0,g,50)
    
    X, Y = np.meshgrid(z, r)
    for z,r in zip(np.ravel(X), np.ravel(Y)):
        print z*1000,r*1000,SpaceChargeField(r,0,z,0,0,0.5*mm)*0.001



#Plot_field()
__name__ = 'sce'
#Nbin = 500
#Ein = 50000 #V/cm
#DiffT = 0.011971
#g = 0.2 #cm
#Dx = 0.2/Nbin
#
#grid = pl.zeros(Nbin)
#Etot = [Ein for i in range(1,Nbin+1)]
#Esc = []
#for i in range(500/2-100,500/2+100):
#    grid[i] = np.random.randint(0,5e5)
#
#    
#
#r = pl.linspace(-0.2,0.2,50)
#z = pl.linspace(0,0.2,50)
#
#for R in r:
#    for Z in z:
#        print Z,R,SpaceChargeField(R,0,Z,0,0,0.1)
#
#
#
#
##print grid
#
#
##parcours de la grid (z)
#for z in range(Nbin):
#    field = 0
#    for zp in range(1,Nbin):
#        sigma = DiffT*np.sqrt(zp*Dx)
#        #print sigma
#        rp = np.random.normal(0,sigma)
#        #print rp
#        if grid[zp] > 0:
#            field += -grid[zp] * SpaceChargeField(0,0,z*Dx,rp,0,zp*Dx)
#    Esc.append(field)
#
#
#Ef = []
#for i in range(Nbin):
#    print Esc[i], Etot[i] + Esc[i], grid[i]
#    Ef.append(Etot[i]+Esc[i])
#plt.plot(Ef)
#plt.show()