
#pragma once

#include <iostream>
#include <cmath>

/*-------------------------------------\
|	integration formulas               |
\-------------------------------------*/

// the integration routine
template<typename Method, typename F, typename Float>
double integrate(F f, void* par, Float a, Float b, int steps, Method m){
	double s = 0;
	double h = (b-a)/steps;
	for (int i = 0; i < steps; ++i)
		s += m(f, par, a + h*i, h);
	return h*s;
}

template<typename F, typename Float>
class iu_transform{
	public:
		iu_transform (F f, Float a) : _f( f ), _a( a ) {}
	
		double operator()(double t, void* par) const{
			double x = _a + (1-t)/t;
			return _f(x,par)/(t*t);
		}
	private:
		F _f;
		Float _a;
};

template<typename Method, typename F, typename Float>
double integrate_iu(F f, void* par, Float a, int steps, Method m){
	iu_transform<F, Float> iu_fun( f, a );
	
	return integrate(iu_fun, par, 0., 1., steps, m);
}
 
// methods
class rectangular{
	public:
	enum position_type { left, middle, right };
	rectangular(position_type pos): position(pos) {}
	template<typename F, typename Float>
	double operator()(F f, void* par, Float x, Float h) const{
		switch(position){
			case left:
				return f(x,par);
			case middle:
				return f(x+h/2,par);
			case right:
				return f(x+h,par);
			}
	}
	private:
	const position_type position;
};
 
class trapezium{
	public:
	template<typename F, typename Float>
	double operator()(F f, void* par, Float x, Float h) const {
		return (f(x,par) + f(x+h,par))/2;
	}
};
 
class simpson{
	public:
	template<typename F, typename Float>
	double operator()(F f, void* par, Float x, Float h) const {
		//std::cout << (f(x,par) + 4*f(x+h/2,par) + f(x+h,par))/6 << " " << x << std::endl;
		return (f(x,par) + 4*f(x+h/2,par) + f(x+h,par))/6;
	}
};

