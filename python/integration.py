import numpy as np
import scipy as sp
from scipy.integrate import quad
import scipy.constants
#from math import np.exp, cos, sqrt

mm = 1.e-3
cm = 1.e-2
g = 0.
q = 0.
p = 0.
eps0 = scipy.constants.epsilon_0
eps1 = 0
eps2 = 0.
eps3 = 0.
e0 = scipy.constants.e
fDiffT = 0.

def RadialChargeDistribution(r, l):
	return ( 1./(fDiffT*fDiffT * l) ) * np.exp( -(r*r)/(2*fDiffT*fDiffT * l) )


def Ebar(x, z, l ,zp):
	#x == rp
	#print x
	E = Field(0.,0.,z,x*cm,0.,zp)
	radDist = RadialChargeDistribution(x,l)
	#print x, E, radDist
	res =   radDist * E*0.01 * x  #RadialDistrib en cm-2, rp*drp en cm2, SCField en V/m. rp doit etre en m dans les params de SCField
	return res;

def compute_Ebar(z, l, zp, diffT, u_eps1, u_eps2, u_eps3, u_p, u_q, u_g):
    global fDiffT, eps1, eps2, eps3, p, q, g
    fDiffT = diffT
    eps1 = u_eps1
    eps2 = u_eps2
    eps3 = u_eps3
    p = u_p
    q = u_q
    g = u_g

    res, error = quad(Ebar,0.,sp.Inf,args=(z,l,zp))

    return res

def Pot(r,phi,z,rp,phip,zp):
    P = np.sqrt(r*r - 2*r*rp*np.cos(phi-phip) + rp*rp)
    #integral = Integral(P,z,zp)
    integral = quad(integrand,0,sp.Inf,args=(P,z,zp))[0]
    #print integral
    Q=e0
    pot = ( Q/(4*np.pi*eps2) ) * ( (1./(np.sqrt(P*P+(z-zp)**2))) - ((eps1-eps2)/((eps1+eps2)*np.sqrt(P*P+(z+zp)**2))) - ((eps3-eps2)/((eps3+eps2)*np.sqrt(P*P+(2*g-z-zp)**2))) + \
                                  (1./((eps1+eps2)*(eps2+eps3))) * integral )
    return pot

def Field(r,phi,z,rp,phip,zp):
    h = 0.01*mm
    m1 = ( Pot(r,phi,z+h,rp,phip,zp) - Pot(r,phi,z-h,rp,phip,zp) ) / (2*h)
    m2 = ( Pot(r,phi,z+2*h,rp,phip,zp) - Pot(r,phi,z-2*h,rp,phip,zp) ) / (4*h)
    m3 = ( Pot(r,phi,z+3*h,rp,phip,zp) - Pot(r,phi,z-3*h,rp,phip,zp) ) / (6*h)

    return - (1.5*m1 - (3./5.)*m2 + 0.1*m3)

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

def integrand(k, P, z, zp):
    return sp.special.j0(k*P) * ( R(k,z,zp)/D(k) )

def compute_pot_correction_term(r, phi, z, rp, phip, zp, u_eps1, u_eps2, u_eps3, u_p, u_q, u_g):
    """ Compute the integral of correction term for potential of point charge in a 3-layers geometry """

    P = np.sqrt(r*r - 2*r*rp*np.cos(phi-phip) + rp*rp)

    global eps1, eps2, eps3, p, q, g
    eps1 = u_eps1
    eps2 = u_eps2
    eps3 = u_eps3
    p = u_p
    q = u_q
    g = u_g

    corr_term, error = quad(integrand,0,sp.Inf,args=(P, z, zp))

    return corr_term

def test(text):
    print "successfully called test() with argument "+text
    x = compute_pot_correction_term(0.,0.,1.e-3,0.,0.,1.e-3, 10.*eps0,eps0,10.*eps0,4.*mm,2.*mm,2.*mm)
    return x

#cm
#z = 0.1*cm
#zp = 0.1*cm
#l = 0.15
#print compute_Ebar(z,l,zp,0.011971,10.*eps0,eps0,10.*eps0,4.*mm,2.*mm,2.*mm)
# z = np.linspace(0,2.*mm,100)
# zp = [0.1e-3,0.5e-3,1.e-3]
# I1 = []
# I2 = []
# I3 = []
# I = [I1,I2,I3]
#
# for i in range(len(zp)):
#     for val in z:
#         I[i].append( ( e0/(4*np.pi*10.*eps0) ) * (1./((10.*eps0+10.*eps0)*(10.*eps0+10.*eps0))) * compute_pot_correction_term(0.,0.,val,0.,0.,zp[i], 10.*eps0,eps0,10.*eps0,4.*mm,2.*mm,2.*mm)[0] )
# plt.plot(z,I[0])
# plt.plot(z,I[1])
# plt.plot(z,I[2])
# plt.show()
