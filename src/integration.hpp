
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

// Does not work apparently ...
// ... No, it definetly does not work !
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


/*	From http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#C.2B.2B 
 *	Modified for semi-improper integration [0,+infty] */
namespace Rosetta {
 
    /* Implementation of Gauss-Legendre quadrature
    *  http://en.wikipedia.org/wiki/Gaussian_quadrature
    *  http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
    */
	template <int N>
    class GaussLegendreQuadrature {
    public:
		enum {eDEGREE = N};
 
        /* 	Compute the integral of a functor
        *
        *   @param a    lower limit of integration
        *   @param b    upper limit of integration
        *   @param f    the function to integrate
        *   @param err  callback in case of problems
        */
        template <typename Function>
        double integrate(double a, double b, Function f) {
            double p = (b - a) / 2;
            double q = (b + a) / 2;
            const LegendrePolynomial& legpoly = s_LegendrePolynomial;
 
            double sum = 0;
            for (int i = 1; i <= eDEGREE; ++i) {
                sum += legpoly.weight(i) * f(p * legpoly.root(i) + q);
            }
 
            return p * sum;
        }
        
        /* Modified function integrate for semi-improper interval of the kind [a,+infty] */
		template <typename Function>
        double integrate_iu(double a, Function f, void* funParams) {
			double b = 1.;
            double p = (b - a) / 2;
            double q = (b + a) / 2;
            const LegendrePolynomial& legpoly = s_LegendrePolynomial;
 
            double sum = 0;
            double t = 0;
            for (int i = 1; i <= eDEGREE; ++i) {
				t = p * legpoly.root(i) + q;
                sum += legpoly.weight(i) * f( a + t/(1-t),funParams )/((1-t)*(1-t));
            }
 
            return p * sum;
        }
 
        /* Print out roots and weights for information
        */
        void print_roots_and_weights(std::ostream& out) const {
            const LegendrePolynomial& legpoly = s_LegendrePolynomial;
            out << "Roots:  ";
            for (int i = 0; i <= eDEGREE; ++i) {
                out << ' ' << legpoly.root(i);
            }
            out << '\n';
            out << "Weights:";
            for (int i = 0; i <= eDEGREE; ++i) {
                out << ' ' << legpoly.weight(i);
            }
            out << '\n';
        }
	private:
        /* 	Implementation of the Legendre polynomials that form
        *   the basis of this quadrature
        */
        class LegendrePolynomial {
        public:
            LegendrePolynomial () {
                // Solve roots and weights
                for (int i = 0; i <= eDEGREE; ++i) {
                    double dr = 1;
 
                    // Find zero
                    Evaluation eval(cos(M_PI * (i - 0.25) / (eDEGREE + 0.5)));
                    do {
                        dr = eval.v() / eval.d();
                        eval.evaluate(eval.x() - dr);
                    } while (fabs (dr) > 2e-16);
 
                    this->_r[i] = eval.x();
                    this->_w[i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
                }
            }
 
            double root(int i) const { return this->_r[i]; }
            double weight(int i) const { return this->_w[i]; }
        private:
            double _r[eDEGREE + 1];
            double _w[eDEGREE + 1];
 
            /* 	Evaluate the value *and* derivative of the
            *   Legendre polynomial
            */
            class Evaluation {
            public:
                explicit Evaluation (double x) : _x(x), _v(1), _d(0) {
                    this->evaluate(x);
                }
 
                void evaluate(double x) {
                    this->_x = x;
 
                    double vsub1 = x;
                    double vsub2 = 1;
                    double f     = 1 / (x * x - 1);
 
                    for (int i = 2; i <= eDEGREE; ++i) {
                        this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
                        this->_d = i * f * (x * this->_v - vsub1);
 
                        vsub2 = vsub1;
                        vsub1 = this->_v;
                    }
                }
 
                double v() const { return this->_v; }
                double d() const { return this->_d; }
                double x() const { return this->_x; }
 
				private:
                double _x;
                double _v;
                double _d;
            };
        };
 
        /* 	Pre-compute the weights and abscissae of the Legendre polynomials
        */
        static LegendrePolynomial s_LegendrePolynomial;
    };
 
    template <int N>
    typename GaussLegendreQuadrature<N>::LegendrePolynomial GaussLegendreQuadrature<N>::s_LegendrePolynomial;
}

