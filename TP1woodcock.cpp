#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "utils.hpp"

using namespace std;
const float int_exp_int_sigma = 1.4163275798634898; //1-exp(-.3)+1./3*(exp(-1.5)-exp(-2.7))+exp(.1)

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

float f_gene(float mu)
{ 
    return mu * inv_int_sigma( -log(ranf()/int_exp_int_sigma) ) + source();
}

vector<float> trajs(float mu, int n)
{
    vector<float> parts(n); // n est la taille du tableau parts
    //int i;
    for(int i = 0 ; i < n ; i++)
    {
	parts.at(i) = f_gene(mu);
    }
    return parts;
}

distrib_t denom(const vector<float>& parts, int nb_segs)
{
    distrib_t distrib(nb_segs);
    for (auto const & val : parts) // sans indice
    {
	if (val>=1) distrib.above++; 
	else if (val< 0) distrib.below++;
	else distrib.segs.at( floor(nb_segs*val) )++;
    }
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


    // vector<float> Xtest = linspace(0, 5, 1000);
    // vector<float> Ytest (1000);
    // for (int i=0; i<1000; i++) Ytest[i] = mu * inv_int_sigma( -log(Xtest[i]/int_exp_int_sigma) );
    // plot(Xtest, Ytest);
    
    vector<float> parts = trajs(mu, nb_parts);
    vector<float> X = linspace(0, 1, nb_segs);
    distrib_t distrib = denom(parts, nb_segs);
    plot(X, distrib.segs, "w l");
    cout << "above / below : " << distrib.above << " / " << distrib.below << endl;
    return EXIT_SUCCESS;
}
