#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "utils.hpp"
#include "vector_v2.hpp"

using namespace std;

struct constants_t {
    constants_t() {}
    constants_t(float mu) :
	C1 (1), 
	C2 (1-2/3.*exp(-0.3/mu)), 
	C3 (C2 +2/3.*exp(-1.5/mu)), 
	l (1/C3), 
	y1 (l*(1-exp(-0.3/mu))), 
	y2 (1-l*exp(-1.5/mu))
    {}
    float C1, C2, C3,
	l, y1, y2;
};

float source() { return 0; }

float sigma(float x)
{
    if (x>0.3 && x<0.7) return 3;
    else return 1;
}

float int_sigma(float x)
{
    if (x<0.3) return x;
    else if (x<0.7) return 3*(x-0.3)+0.3;
    else return (x-.7)+1.5;
}

float inv_int_sigma(float y)
{
    if (y<0.3) return y;
    else if (y<1.5) return (y-0.3)/3+0.3;
    else return (y-1.5)+0.7;
}

float P(float x, float mu, const constants_t& cst)
{
    if (x<0.3) return cst.l*exp(-x/mu);
    else if (x<0.7) return cst.l*exp((0.6-3*x)/mu);
    else return cst.l*exp(-(x+0.8)/mu);
}

float F_repartition(float x, float mu, const constants_t& cst)
{
    if (x<0.3) return cst.l*(cst.C1-exp(-x/ mu));
    else if (x<0.7) return cst.l*(cst.C2-1/3.*exp((-3*x+0.6)/mu));
    else return cst.l*(cst.C3-exp(-(x+0.8)/mu));
}

float inv_F_repartition(float y, float mu, const constants_t& cst)
{
    if (y<cst.y1) return -mu * log(cst.C1-y/cst.l);
    else if (y<cst.y2) return 1/3.*( 0.6-mu*log(3*(cst.C2-y/cst.l)) );
    else return -0.8 - mu*log(cst.C3-y/cst.l);
}
    

float f_gene(float mu, const constants_t &cst)
{
    return inv_F_repartition(ranf(), mu, cst) + source();
}

vectorV2<float> trajs(float mu, int n, const constants_t& cst)
{
    vectorV2<float> parts(n); // n est la taille du tableau parts
    //int i;
    for(int i = 0 ; i < n ; i++)
    {
	parts.at(i) = f_gene(mu, cst);
    }
    return parts;
}

distrib_t denom(const vectorV2<float>& parts, int nb_segs)
{
    distrib_t distrib(nb_segs);
    for (auto const & val : parts) // sans indice
    {
	if (val>=1) distrib.above++; 
	else if (val< 0) distrib.below++;
	else distrib.segs.at( floor(nb_segs*val) )++;
    }
    distrib.above = distrib.above * nb_segs / parts.size();
    for (auto &s : distrib.segs) s = s * nb_segs / parts.size();
    distrib.below = distrib.below * nb_segs / parts.size();
    return distrib;
}

int main(int argc, char**argv)
{
    if (argc!=4) return EXIT_FAILURE;
    float mu;
    int nb_parts,
	nb_segs;
    mu = atof(argv[1]);
    nb_parts = atoi(argv[2]); // i pour pas float bro
    nb_segs  = atoi(argv[3]);
    constants_t cst(mu);

    // vector<float> Xtest = linspace(0, 1, 1000);
    // vector<float> Ytest (1000);
    // for (int i=0; i<1000; i++) Ytest[i] = inv_F_repartition(Xtest[i], mu);
    // plot(Xtest, Ytest);

    vectorV2<float> parts = trajs(mu, nb_parts, cst);
    vectorV2<float> X = linspace(0, 1, nb_segs);
    distrib_t distrib = denom(parts, nb_segs);
    plot(X, distrib.segs, "w l");

    vectorV2<float> Py (nb_segs) ;
    for (int i = 0; i<nb_segs; i++) Py.at(i) = P(X.at(i), mu, cst);
    plot(X, Py, "w l");

    cout << "above / below : " << distrib.above << " / " << distrib.below << endl
	 << "diff normalized : " << (distrib.segs-Py).norm()/Py.norm() << endl;
    return EXIT_SUCCESS;
}
