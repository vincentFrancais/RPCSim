import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import scipy as sp
from scipy.integrate import quad
import scipy.constants

mm = 1.e-3
g = q = 2.*mm
p = 4.*mm
eps0 = scipy.constants.epsilon_0
eps1 = eps3 = 10.*eps0
eps2 = eps0
e0 = scipy.constants.e

def D(k):
    return (eps1+eps2)*(eps2+eps3)*(1.-np.exp(-2.*k*(p+q))) \
            - (eps1-eps2)*(eps2+eps3)*(np.exp(-2.*k*p)-np.exp(-2.*k*q)) \
            - (eps1+eps2)*(eps2-eps3)*(np.exp(-2.*k*(p-g))-np.exp(-2.*k*(q+g))) \
            + (eps1-eps2)*(eps2-eps3)*(np.exp(-2.*k*g)-np.exp(-2.*k*(p+q-g)))

def R(k,z,zp):
    return (eps1+eps2)**2 * (eps2+eps3)**2 * (np.exp(k*(-2.*p-2.*q+z-zp)) + np.exp(k*(-2.*p-2.*q-z+zp))) \
            - (eps1+eps2)**2 * (eps2-eps3)**2 * np.exp(k*(-4.*g-2.*q+z+zp)) \
            - 4.*eps1*eps2*(eps2+eps3)**2 * np.exp(k*(-2.*q-z-zp)) - (eps1-eps2)**2 * (eps2+eps3)**2 * np.exp(k*(-2.*p-z-zp)) \
            - (eps1**2 - eps2**2)*(eps2-eps3)**2 * np.exp(k*(-4.*g+z+zp)) \
            + (eps1**2 - eps2**2)*(eps2+eps3)**2 * ( -np.exp(k*(-2.*p-2.*q-z-zp)) + np.exp(k*(-2.*p+z-zp)) + np.exp(k*(-2.*p-z+zp)) ) \
            - 4.*(eps1**2 - eps2**2)*eps2*eps3*np.exp(k*(-2.*p-2.*q+z+zp)) - 4.*(eps1+eps2)**2 * eps2*eps3*np.exp(k*(-2.*p+z+zp)) \
            + (eps1-eps2)**2 * (eps2**2 - eps3**2) * np.exp(k*(-2.*g-z-zp)) + 4.*eps1*eps2*(eps2**2 - eps3**2)*np.exp(k*(2.*g-2.*p-2.*q-z-zp)) \
            + (eps1+eps2)**2 * (eps2**2 - eps3**2) * ( -np.exp(k*(-2.*g-2.*q+z-zp)) - np.exp(k*(-2.*g-2.*q-z+zp)) + np.exp(k*(-2.*g-2.*p-2.*q+z+zp)) ) \
            + (eps1**2 - eps2**2)*(eps2**2 - eps3**2) * ( np.exp(k*(-2.*g-2.*q-z-zp)) - np.exp(k*(-2.*g+z-zp)) - np.exp(k*(-2.*g-z+zp)) + np.exp(k*(-2.*g-2.*p+z+zp)) )

def integrand(k,P,z,zp):
    return sp.special.j0(k*P) * ( R(k,z,zp)/D(k) )

def Integral(P,z,zp):
    I = 0.
    h = 1.1
    n = 0
    epsilon = 1e-8
    while(True):
        a = n*h
        b = (n+1)*h
        x = (a + b)*0.5
        previous = I
        val = ((b - a)/6) * (integrand(a,P,z,zp)+4*integrand(x,P,z,zp)+integrand(b,P,z,zp))
        #print val,a,b
        if abs(val) == sp.Inf:
            n+=1
            continue
        I += val

        #print I,previous, abs(I - previous)
        if n>300:# and abs(I - previous) < epsilon:
            break

        n+=1

    return I


def Pot(r,phi,z,rp,phip,zp):
    P = np.sqrt(r*r - 2*r*rp*np.cos(phi-phip) + rp*rp)
    #integral = Integral(P,z,zp)
    integral = quad(integrand,0,sp.Inf,args=(P,z,zp))[0]
    #print integral
    Q=e0
    pot = ( Q/(4*np.pi*eps2) ) * ( (1./(np.sqrt(P*P+(z-zp)**2))) - ((eps1-eps2)/((eps1+eps2)*np.sqrt(P*P+(z+zp)**2))) - ((eps3-eps2)/((eps3+eps2)*np.sqrt(P*P+(2*g-z-zp)**2))) + \
                                  (1./((eps1+eps2)*(eps2+eps3))) * integral )
    return pot

def Plot_Pot():
    r = pl.linspace(-1.*mm,1.*mm,40)
    z = pl.linspace(0,g,40)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.gca().patch.set_facecolor('white')
    ax.w_xaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_yaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_zaxis.set_pane_color((1., 1., 1., 1.0))
    #x = y = np.arange(-3.0, 3.0, 0.05)
    X, Y = np.meshgrid(z, r)
    zs = np.array([Pot(r,0.,z,0.,0.,0.5*mm) for z,r in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)

    ax.set_xlabel(r'$z\,[mm]$', fontsize=20)
    ax.set_ylabel(r'$r\,[mm]$', fontsize=20)
    ax.set_zlabel(r'$\Phi\,[10^{-5}\,V]$', fontsize=20)

    ax.plot_wireframe(X*1000, Y*1000, Z*1e5)
    #ax.set_zlim(0,1.e-8)
    #ax.dist=9
    plt.savefig("pot.pdf")
    plt.show()

def Field(r,phi,z,rp,phip,zp):
    h = 0.01*mm
    m1 = ( Pot(r,phi,z+h,rp,phip,zp) - Pot(r,phi,z-h,rp,phip,zp) ) / (2*h)
    m2 = ( Pot(r,phi,z+2*h,rp,phip,zp) - Pot(r,phi,z-2*h,rp,phip,zp) ) / (4*h)
    m3 = ( Pot(r,phi,z+3*h,rp,phip,zp) - Pot(r,phi,z-3*h,rp,phip,zp) ) / (6*h)

    return - (1.5*m1 - (3./5.)*m2 + 0.1*m3)
    #return ( Pot(r,phi,z+h,rp,phip,zp) - Pot(r,phi,z-h,rp,phip,zp) ) / 2*h

def plot_Field():
    r = pl.linspace(-1.*mm,1.*mm,40)
    z = pl.linspace(0,g,40)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.gca().patch.set_facecolor('white')
    ax.w_xaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_yaxis.set_pane_color((1., 1., 1., 1.0))
    ax.w_zaxis.set_pane_color((1., 1., 1., 1.0))
    #x = y = np.arange(-3.0, 3.0, 0.05)
    X, Y = np.meshgrid(z, r)
    zs = np.array([Field(r,0.,z,0.,0.,0.5*mm) for z,r in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)

    ax.set_xlabel(r'$z\,[mm]$', fontsize=20)
    ax.set_ylabel(r'$r\,[mm]$', fontsize=20)
    ax.set_zlabel(r'$E\,[V/m]$', fontsize=20)

    ax.plot_wireframe(X*1000, Y*1000, Z)
    #ax.set_zlim(0,1.e-8)
    #ax.dist=9
    plt.savefig("field.pdf")
    plt.show()

def Plot_integrand():
    z = 1.* mm
    zp = 0.48 * mm
    P = 0.
    k = pl.linspace(0.0001,3.,10000)
    I =  integrand(k*1.e3,P,z,zp)
    r = R(k,z,zp)
    d = D(k)
    plt.plot(k,I)
    plt.xlabel(r'$\kappa$', fontsize=20)
    plt.ylabel(r'$J_0(P\kappa) \, R(\kappa,z,z\textasciiacute  )/D(\kappa)$', fontsize=20)
    plt.grid(True)
    plt.savefig('plot_integrand.pdf')
    plt.show()

def Plot_integral():
    zp = [0.1e-3,0.5e-3,1.e-3]
    z = pl.linspace(0,g,100)
    P = 0.
    I1 = []
    I2 = []
    I3 = []
    I = [I1,I2,I3]
    print quad(integrand,0.,sp.Inf,args=(P,1.e-3,0.48e-3))[0]
    print quad(integrand,0.,3.,args=(P,1.e-3,0.48e-3))[0]
    for i in range(len(zp)):
        for val in z:
            #x = Integral(P,val,zp)
            x = (e0/(4*np.pi*eps2) )  *  (1./((eps1+eps2)*(eps2+eps3))) * quad(integrand,0,sp.Inf,args=(P,val,zp[i]))[0]
        #print x
            I[i].append(x)
    plt.plot(z*1e3,I[0],label='z=0.1')
    plt.plot(z*1e3,I[1],'r',label='z=0.5')
    plt.plot(z*1e3,I[2],label='z=1.')
    plt.xlabel(r'$z[m]$', fontsize=20)
    plt.grid(True)
    plt.legend()
    #plt.savefig('plot_integral_pot.pdf')
    plt.show()

def compute_terms(r,phi,z,rp,phip,zp):
    P = np.sqrt(r*r - 2*r*rp*np.cos(phi-phip) + rp*rp)
    term1 = (1./(np.sqrt(P*P+(z-zp)**2)))
    term2 = ((eps1-eps2)/((eps1+eps2)*np.sqrt(P*P+(z+zp)**2)))
    term3 = ((eps3-eps2)/((eps3+eps2)*np.sqrt(P*P+(2*g-z-zp)**2)))
    #term4 = (1./((eps1+eps2)*(eps2+eps3))) * integral

    return term1, term2, term3 #term4

def Plot_terms():
    r = 0. * mm
    rp = 0. * mm
    phi = phip = 0.
    #zp = [0.1e-3,0.5e-3,1.e-3]
    zp = 1. * mm
    z = pl.linspace(0.,g,1000)
    t = []
    for val in z:
        term1, term2, term3 = compute_terms(r,phi,val,rp,phip,zp)
        t.append( ( e0/(4*np.pi*eps2) ) * term1)

    t = np.array(t)
    plt.plot(z,t)
    #plt.ylim(0.,5)
    plt.show()

#Plot_integrand()
#Plot_Pot()
Plot_integral()
#Plot_terms()

#plot_Field()
