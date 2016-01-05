import numpy as np


def legendre_poly(n,x):
	if n==0:
		return x*0. + 1.
	
	elif n==1:
		return x
	
	else:
		return ( (2.*n - 1) * x * legendre_poly(n-1,x) - (n-1) * legendre_poly(n-2,x) ) / n
		
def legendre_poly_der(n,x):
	if n==0:
		return x*0.
	
	elif n==1:
		return 1.
	
	else:
		return ( n/(x**2 - 1) ) * ( x*legendre_poly(n,x) - legendre_poly(n-1,x) )
		
def legendre_roots(order,tolerance=1e-20):
	roots = []
	err = 0
	if order<2:
		err = 1
		print 'error leg roots'
		return [roots,err]

	for i in range(1,int(order)/2 +1):
		x = np.cos( np.pi * (i-0.25)/(order+0.5) )
		error = 0.
		iters = 0
		
		while (error>tolerance):
			dx = legendre_poly(order,x)/legendre_poly_der(order,x)
			x -= dx
			iters += 1
			error = dx
			if(iters>10000):
				break
		roots.append(x)
		
	roots = np.array(roots)
	if order%2==0:
		roots = np.concatenate( (-1.0*roots, roots[::-1]) )
	else:
		roots = np.concatenate( (-1.0*roots, [0.0], roots[::-1]) )
		
	err = 0
	return [roots,err]
	
def gaussLegendre_weights(order):
	w = []
	xi,err = legendre_roots(order)
	
	if err!=0:
		return [[],[],err]
	
	w = 2. / ( (1. - xi**2)*legendre_poly_der(order,xi)**2 )
	
	return [w,xi,0]
	

#w,xi,err = gaussLegendre_weights(30)
#print w
#print xi
#print err

from numpy.polynomial.laguerre import *
#186
out = open("gauss-laguerre.hpp", "w")
#for i in range(89,90):
i = 120
print i
x,w = laggauss(i)
out.write("double x"+str(i)+"["+str(i)+"] = {")
for j in range(len(x)-1):
	out.write( str(x[j])+"," )
out.write( str(x[len(x)-1])+"};\n" )

out.write("double w"+str(i)+"["+str(i)+"] = {")
for j in range(len(w)-1):
	out.write( str(w[j])+"," )
out.write( str(w[len(w)-1])+"};\n\n" )
