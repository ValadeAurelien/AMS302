#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
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

// Structure pour garder en mémoire le départ à l'arrivée de la particule
struct part_trajectory_t
{
    part_trajectory_t () {};
    part_trajectory_t (float _b, float _e) : begin(_b), end(_e) {};

    float begin,
	end;
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

float Phi(float x, float mu, const constants_t& cst)
{
    // if (x<0.3) return cst.l/mu*exp(-x/mu);
    // else if (x<0.7) return cst.l/mu*exp((0.6-3*x)/mu);
    // else return cst.l/mu*exp(-(x+0.8)/mu);
    if (x<0.3) return 1/mu*exp(-x/mu);
    else if (x<0.7) return 1/mu*exp((0.6-3*x)/mu);
    else return 1/mu*exp(-(x+0.8)/mu);
}

float F_repartition_propagateur(float x, float mu, const constants_t& cst)
{
    if (x<0.3) return cst.l*(cst.C1-exp(-x/ mu));
    else if (x<0.7) return cst.l*(cst.C2-1/3.*exp((-3*x+0.6)/mu));
    else return cst.l*(cst.C3-exp(-(x+0.8)/mu));
}

float inv_F_repartition_propagateur(float y, float mu, const constants_t& cst)
{
    if (y<cst.y1) return -mu * log(cst.C1-y/cst.l);
    else if (y<cst.y2) return 1/3.*( 0.6-mu*log(3*(cst.C2-y/cst.l)) );
    else return -0.8 - mu*log(cst.C3-y/cst.l);
}
    
/* Variable aléatoire suivant la densité de probabilité
   du libre parcours d'un neutron */
part_trajectory_t one_traj(float mu, const constants_t& cst)
{
    float src = source();
    return part_trajectory_t( src, inv_F_repartition_propagateur(ranf(), mu, cst) + src);
}

vector<part_trajectory_t> trajs(float mu, int n, const constants_t& cst)
{
    vector<part_trajectory_t> parts(n); // n est la taille du tableau parts
    //int i;
    for(int i = 0 ; i < n ; i++)
    {
	parts.at(i) = one_traj(mu, cst);
    }
    return parts;
}

distrib_t denom(const vector<part_trajectory_t>& parts, int nb_segs, const constants_t& cst)
{
    distrib_t distrib(nb_segs);
    for (auto const & val : parts) // sans indice
    {
	if (val.end>=1) distrib.above++; 
	else if (val.end< 0) distrib.below++;
	else distrib.segs.at( floor(nb_segs*val.end) )++;
    }
    distrib.above = distrib.above * nb_segs / (cst.l*parts.size()); // on divise par l'intégrale du propagateur
    distrib.segs  = distrib.segs  * nb_segs / (cst.l*parts.size());
    distrib.below = distrib.below * nb_segs / (cst.l*parts.size());
    return distrib;
}

int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc<5) {
	cout << "Enter mu, nb_particules, nb_segments, "
	     << "output style (none=1, plot=2, file=3, all=4) [, filename]"
	     << endl;
	return EXIT_FAILURE;
    }
    float mu = atof(argv[1]);                 // mu du problème
    int nb_parts = floor(atof(argv[2])),    // Nombre de particules de l'échantillon (floor(atof()) pour pouvoir utiliser 1e6 etc...
	nb_segs  = floor(atof(argv[3]));    // Finesse de segmentation de l'intervalle pour calculer le flux
    int output_style = atoi(argv[4]);
    if ( !mu ||
	 !nb_parts ||
	 !nb_segs ||
	 !output_style ) {
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }
    string fname;
    if (output_style>2) {
	if (argc>5)
	    fname = argv[5];
	else
	    fname = "output_TP1_sigmavar";
    }

    // constants du problème
    constants_t cst(mu);
    // Courbe théorique
    vectorV2<float> X = linspace(0, 1, nb_segs),
	Py (nb_segs) ;
    for (int i = 0; i<nb_segs; i++)
	Py.at(i) = Phi(X.at(i), mu, cst);
    
    // Monte Carlo
    vector<part_trajectory_t> parts = trajs(mu, nb_parts, cst);
    distrib_t distrib = denom(parts, nb_segs, cst);

    // for (int i=0; i<nb_segs-1; i++)
    // 	cout << nb_segs*(distrib.segs.at(i+1)-distrib.segs.at(i)) << endl;
    // Affichage
    if (output_style>2) {
	ofstream file (fname, fstream::out);
	file << "#";
	for (int i=0; i<argc; i++) file << argv[i] << " ";
	file << endl;
	file << "#distrib" << endl << "#X     Phi(MC)     Phi(th)" << endl;
	for (int i=0; i<nb_segs; i++)
	    file << X.at(i) << " " << distrib.segs.at(i) << " " << Py.at(i) << endl;
	file << "#(above 1) / (below 0) : "
	     << distrib.above << " / " << distrib.below << endl
	     << "#diff normalized : " << (distrib.segs-Py).norm()/Py.norm()
	     << endl;
	file.close();
    }
    if (!(output_style % 2)) {
	plot(X, distrib.segs, "w l", "set title 'MC'; set yrange [0:]; ");
	plot(X, Py, "w l", "set title 'théorique'; set yrange [0:]; ");
	cout << "#(above 1) / (below 0) : " << distrib.above << " / " << distrib.below << endl
	 << "#diff normalized : " << (distrib.segs-Py).norm()/Py.norm() << endl;
    }
    return EXIT_SUCCESS;
}
