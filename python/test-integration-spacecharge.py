import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy as sp
from scipy.integrate import quad


g = q = 2.
p = 4.
eps0 = 8.854e-12
eps1 = eps3 = 10.*eps0
eps2 = eps0
e0 = 1.60217657e-19

def D(k):
    return (eps1+eps2)*(eps2+eps3)*(1-np.exp(-2*k*(p+q))) - (eps1-eps2)*(eps2+eps3)*(np.exp(-2*k*p)-np.exp(-2*k*q)) - (eps1+eps2)*(eps2-eps3)*(np.exp(-2*k*(p-g))-np.exp(-2*k*(q+g))) + (eps1-eps2)*(eps2-eps3)*(np.exp(-2*k*g)-np.exp(-2*k*(p+q-g)))

def R(k,z,zp):
    return (eps1+eps2)**2 * (eps2+eps3)**2 * (np.exp(-2*p-2*q+z-zp) + np.exp(-2*p-2*q-z+zp)) \
            - (eps1+eps2)**2 * (eps2-eps3)**2 * np.exp(k*(-4*g-2*q+z+zp)) \
            + 4*eps1*eps2*(eps2+eps3)**2 * np.exp(k*(-2*q-z-zp)) - (eps1-eps2)**2 * (eps2+eps3)**2 * np.exp(-k*(-2*p-z-zp)) \
            - (eps1**2 - eps2**2)*(eps2-eps3)**2 * np.exp(k*(-4*g+z+zp)) \
            + (eps1**2 - eps2**2)*(eps2+eps3)**2 * ( -np.exp(k*(-2*p-2*q-z-zp)) + np.exp(k*(-2*p+z-zp)) + np.exp(k*(-2*p-z+zp)) ) \
            - 4*(eps1**2 - eps2**2)*eps2*eps3*np.exp(k*(-2*p-2*q+z+zp)) - 4*(eps1**2 + eps2**2)*eps2*eps3*np.exp(k*(-2*p+z+zp)) \
            + (eps1-eps2)**2 * (eps2**2 - eps3**2) * np.exp(k*(-2*g-z-zp)) + 4*eps1*eps2*(eps2**2 - eps3**2)*np.exp(k*(-2*g-2*p-2*q-z-zp)) \
            + (eps1+eps2)**2 * (eps2**2 - eps3**2) * ( -np.exp(k*(-2*g-2*q+z-zp)) - np.exp(k*(-2*g-2*q-z+zp)) + np.exp(k*(-2*g-2*p-2*q+z+zp)) ) \
            + (eps1**2 - eps2**2)*(eps2**2 - eps3**2) * ( np.exp(k*(-2*g-2*q-z-zp)) - np.exp(k*(-2*g+z-zp)) - np.exp(k*(-2*g-z+zp)) + np.exp(k*(-2*g-2*p+z+zp)) )

def integrand(k,P,z,zp):
    return sp.special.j0(k*P) * R(k,z,zp)/D(k)
    
def Pot(r,phi,z,rp,phip,zp):
    P = r*r - 2*r*rp*np.cos(phi-phip) + rp*rp
    integral = 0.#quad(integrand,0,sp.Inf,args=(P,z,zp))[0]
    #print integral
    Q=e0
    pot = ( Q/(4*np.pi*eps2) ) * ( (1./(np.sqrt(P*P+(z-zp)**2))) - ((eps1-eps2)/((eps1+eps2)*np.sqrt(P*P+(z+zp)**2))) - ((eps3-eps2)/((eps3+eps2)*np.sqrt(P*P+(2*g-z-zp)**2))) + \
                                  (1./((eps1+eps2)*(eps2+eps3))) * integral )
    return pot
    
    
r = pl.linspace(-1.,1.,30)
z = pl.linspace(0,g,30)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#x = y = np.arange(-3.0, 3.0, 0.05)
X, Y = np.meshgrid(z, r)
zs = np.array([Pot(r,0.,z,0.,0.,0.5) for z,r in zip(np.ravel(X), np.ravel(Y))])
Z = zs.reshape(X.shape)

ax.set_xlabel('z')
ax.set_ylabel('r')
ax.set_zlabel('E')

ax.plot_wireframe(X, Y, Z)
#ax.set_zlim(0,1.e-8)
#ax.dist=9
plt.savefig("pot.pdf")
plt.show()