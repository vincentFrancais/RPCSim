#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>

#include "gsl/gsl_sf_laguerre.h"

#define M_PI  3.14159265358979323846

using namespace std;

double legendre_poly(int n, double x){
	if (n==0)
		return x*0. + 1.;
		
	else if (n==1)
		return x;
		
	else
		return ( (2.*n - 1.) * x * legendre_poly(n-1,x) - (n-1) * legendre_poly(n-2,x) ) / n;
}

double laguerre_poly(int n, double x){
	if (n==0)
		return x*0. + 1.;
		
	else if (n==1)
		return -1.*x + 1;
		
	else
		return ( (2.*n - 1. - x) * laguerre_poly(n-1,x) - (n-1) * laguerre_poly(n-2,x) ) / n;
}

double legendre_poly_der(int n, double x){
	if (n==0)
		return x*0.;
	
	else if (n==1)
		return 1.;
	
	else
		return ( n/(x*x - 1.) ) * ( x*legendre_poly(n,x) - legendre_poly(n-1,x) );
}

double laguerre_poly_der(int n, double x){
	if (n==0)
		return x*0.;
	
	else if (n==1)
		return -1.;
	
	else
		return ( n/x ) * ( laguerre_poly(n,x) - laguerre_poly(n-1,x) );
}
		
vector<double> legendre_roots(int order){
	double tolerance = 1.e-20;
	vector<double> roots;
	if (order<2)
		return roots;

	for (int i=1; i<int(order*0.5+1); i++){
		double x = cos( M_PI * (i-0.25)/(order+0.5) );
		double error = 0.;
		int iters = 0;
		double dx = 0.;
		double prev;
		while (error>tolerance){
			dx = legendre_poly(order,x)/legendre_poly_der(order,x);
			x -= dx;
			iters++;
			error = dx;
			if (iters>1 and abs(prev-dx)<=tolerance)
				break;
			prev = dx;
			if (iters > 10000)
				break;
		}
		roots.push_back(x);
	}
	
	vector<double> roots_sym(roots);
	//vector<double> roots_final;
	
	for(uint i=0; i<roots.size(); i++)
		roots_sym[i] = -1.*roots[i];
		
	if (order%2==0){
		//roots_final.resize(2*roots.size());
		//merge ( roots_sym.begin(), roots_sym.end(), roots.begin(), roots.end(), roots_final.begin() );
		reverse(roots.begin(), roots.end());
		roots_sym.insert( roots_sym.end(), roots.begin(), roots.end() );
	}
	else{
		//roots_final.resize(2*roots.size()+1);
		roots_sym.resize(roots_sym.size()+1,0.0);
		reverse(roots.begin(), roots.end());
		//merge ( roots_sym.begin(), roots_sym.end(), roots.begin(), roots.end(), roots_final.begin() );
		roots_sym.insert( roots_sym.end(), roots.begin(), roots.end() );
	}
		
	//err = 0
	return roots_sym;
}

vector<double> laguerre_roots(int order){
	double tolerance = 1.e-20;
	vector<double> roots;
	if (order<2)
		return roots;

	for (int i=1; i<order; i++){
		double x = cos( M_PI * (i-0.25)/(order+0.5) );
		double error = 0.;
		int iters = 0;
		double dx = 0.;
		double prev;
		while (error>tolerance){
			dx = laguerre_poly(order,x)/laguerre_poly_der(order,x);
			x -= dx;
			iters++;
			error = dx;
			if (iters>1 and abs(prev-dx)<=tolerance)
				break;
			prev = dx;
			if (iters > 10000)
				break;
		}
		roots.push_back(x);
	}
		
	//err = 0
	return roots;
}

int gaussLegendre_weights(int order, vector<double>& w, vector<double>& xi){
	w.clear();
	xi.clear();
	xi = legendre_roots(order);

	if (xi.size() < 1)
		return 1;
	
	for (vector<double>::iterator it=xi.begin(); it!=xi.end(); ++it)
		w.push_back( 2. / ( (1. - (*it)*(*it))*(legendre_poly_der(order,*it)*legendre_poly_der(order,*it)) ) );
	
	return 0;
}

int gaussLaguerre_weights(int order, vector<double>& w, vector<double>& xi){
	w.clear();
	xi.clear();
	xi = laguerre_roots(order);

	if (xi.size() < 1)
		return 1;
	
	for (vector<double>::iterator it=xi.begin(); it!=xi.end(); ++it)
		w.push_back( *it / ( (order+1)*(order+1) * laguerre_poly(order+1,*it) * laguerre_poly(order+1,*it) ) );
	
	return 0;
}

double lege_inte(double (*f)(double), double a, double b, int order)
{
	double c1 = (b - a) / 2, c2 = (b + a) / 2, sum = 0;
	vector<double> w,x;
	int err = gaussLegendre_weights(order,w,x);
	for (int i = 0; i < order; i++)
		sum += w[i] * f(c1 * x[i] + c2);
	return c1 * sum;
}

double lag_inte(double (*f)(double), int order)
{
	double sum = 0;
	//vector<double> w,x;
	//int err = gaussLaguerre_weights(order,w,x);
	double x[120] = {0.0119983762543,0.0632209692953,0.155383485473,0.28852206346,0.462661866639,0.677833155992,0.934072827013,1.23142472545,1.56993976333,1.94967599191,2.37069866639,2.83308031237,3.33690079718,3.88224740774,4.46921493546,5.09790576866,5.76842999291,6.48090549947,7.2354581022,8.03222166319,8.87133822756,9.75295816755,10.6772403365,11.6443522328,12.654470175,13.7077794871,14.8044746963,15.9447597422,17.128848199,18.3569635108,19.6293392412,20.946219337,22.3078584077,23.7145220214,25.1664870174,26.6640418374,28.2074868756,29.7971348489,31.4333111889,33.1163544568,34.8466167819,36.6244643264,38.4502777774,40.3244528677,42.2474009279,44.2195494722,46.2413428186,48.3132427474,50.4357292011,52.6093010263,54.8344767633,57.1117954857,59.4418176938,61.8251262658,64.2623274721,66.7540520571,69.3009563933,71.9037237156,74.5630654404,77.2797225773,80.0544672432,82.888104284,85.781473018,88.7354491082,91.7509465782,94.8289199836,97.9703667549,101.176329728,104.447899881,107.7862193,111.192484394,114.667949393,118.213930145,121.831808266,125.523035667,129.2891395,133.131727591,137.052494396,141.053227565,145.13581517,149.30225371,153.554656964,157.895265845,162.326459358,166.850766855,171.470881753,176.189676957,181.010222244,185.935803938,190.969947245,196.116441719,201.379370418,206.763143427,212.27253659,217.912736487,223.689392951,229.608680739,235.677372427,241.902925156,248.293584619,254.858510746,261.607930947,268.553328794,275.707678808,283.085742141,290.704443868,298.583361641,306.745369267,315.217500686,324.032135617,333.22866878,342.855931753,352.97583572,363.669096393,375.044737412,387.257015558,400.538537634,415.274303677,432.205390729,453.253394274};
	double w[120] = {0.030424717619,0.0672913283193,0.0964352506624,0.115129252168,0.122552687271,0.119668016561,0.108763809885,0.092794638375,0.0747167210747,0.0569804406536,0.0412606143417,0.0284206711032,0.0186468081948,0.0116650258488,0.00696327443825,0.0039686793019,0.00216063811099,0.00112402024527,0.000558904907562,0.000265680083214,0.000120753056403,5.24801680556e-05,2.18108359466e-05,8.66834132317e-06,3.29442235584e-06,1.19724484283e-06,4.16022653944e-07,1.38210533295e-07,4.38940881341e-08,1.3324580055e-08,3.8656035949e-09,1.07157663133e-09,2.83785100154e-10,7.17839538101e-11,1.73396543152e-11,3.99876220822e-12,8.80182347644e-13,1.84869207066e-13,3.70405422781e-14,7.07750112572e-15,1.28924187679e-15,2.23818294918e-16,3.7017845006e-17,5.83070856439e-18,8.74296909333e-19,1.24752426462e-19,1.69320146558e-20,2.18497677493e-21,2.6795640924e-22,3.12142218738e-23,3.45219122109e-24,3.62297023797e-25,3.60600437631e-26,3.40199150022e-27,3.04038914621e-28,2.57244850094e-29,2.05924014356e-30,1.55854533288e-31,1.11450220774e-32,7.52446056235e-34,4.79261984814e-35,2.87759488857e-36,1.62736954526e-37,8.66098912782e-39,4.33392867444e-40,2.03713352294e-41,8.98571306274e-43,3.71563467312e-44,1.43877223943e-45,5.21121094169e-47,1.76342998742e-48,5.56815205104e-50,1.63844123318e-51,4.48665528589e-53,1.14172441513e-54,2.69580675329e-56,5.89674495608e-58,1.19289391115e-59,2.22785374114e-61,3.83400781992e-63,6.06792249206e-65,8.81325551045e-67,1.1721297997e-68,1.42407268906e-70,1.57657605393e-72,1.58621642276e-74,1.44622450617e-76,1.19127179958e-78,8.83630358921e-81,5.88159078053e-83,3.49984860364e-85,1.85427185619e-87,8.70895021238e-90,3.6088390452e-92,1.31263126425e-94,4.16731507206e-97,1.14774987808e-99,2.72393851427e-102,5.52960211121e-105,9.523254877e-108,1.37888288005e-110,1.66156257539e-113,1.64744464962e-116,1.32686646948e-119,8.55493529604e-123,4.34212400096e-126,1.70164756662e-129,5.03376467215e-133,1.09437279934e-136,1.69353193747e-140,1.79426234779e-144,1.24028901742e-148,5.26267976582e-153,1.26600031832e-157,1.55119641842e-162,8.32117793288e-168,1.55644029687e-173,6.97371948576e-180,3.64809173124e-187,3.52524257384e-196};


	for (int i = 0; i < order; i++)
		sum += w[i] * f(x[i]) * exp(x[i]);
	return sum;
}

double fun(double x){
	//return exp(-x);
	return 1./((x+1.)*sqrt(x));
}

int main(){
	vector<double> w,x;
	int err = gaussLegendre_weights(5,w,x);
	cout << err << endl;
	for (std::vector<double>::iterator it=x.begin(); it!=x.end(); ++it)
		cout << *it << "\t";
	cout << endl;
	for (std::vector<double>::iterator it=w.begin(); it!=w.end(); ++it)
		cout << *it << "\t";
		
	cout << endl;
	cout << lege_inte(exp, -3, 3, 5) << endl;
	cout << exp(3) - exp(-3) << endl;
	
	cout << "---------" <<endl;
	cout << gsl_sf_laguerre_n(5,0,0.85) <<endl;
	cout << laguerre_poly(5,0.85) <<endl;
	cout << "---------" <<endl;
	cout << lag_inte(&fun,120)<<endl;
	
	
	/*
	int err2 = gaussLaguerre_weights(5,w,x);
	for (std::vector<double>::iterator it=x.begin(); it!=x.end(); ++it)
		cout << *it << "\t";
	cout << endl;
	for (std::vector<double>::iterator it=w.begin(); it!=w.end(); ++it)
		cout << *it << "\t";
*/
}
