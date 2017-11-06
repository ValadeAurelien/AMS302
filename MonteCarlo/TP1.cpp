#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "vector_v2.hpp"
#include "utils.hpp"
#include <stdexcept>

using namespace std;

// Deux types de source
enum source_t
{
    CSTE = 0,
    DELTA = 1
}; 

// Structure pour garder en mémoire le départ à l'arrivée de la particule
struct part_trajectory_t
{
    part_trajectory_t () {};
    part_trajectory_t (float _b, float _e) : begin(_b), end(_e) {};

    float begin,
	end;
};

/* Terme source, suivant qu'il est uniforme (typ_source vrai)
                  ou non-nul seulement en 0 (typ_source faux) */
float source(source_t stype)
{
    switch(stype){
    case CSTE :
	return ranf();
    case DELTA:
	return 0;
    }
}

/* Phi théorique */
float Phi(float x, float mu, float sigma, source_t stype)
{
    switch(stype){
    case CSTE :
	if (mu>0)
	    return 1/sigma*(1-exp(-sigma/mu*x));
	else
	    return 1/sigma*(1-exp(-sigma/mu*(x-1)));
    case DELTA :
	return 1/mu*exp(-sigma/mu*x);
    }
}

/* Fonction qui renvoie la distance en fonction de la 
   proba uniforme sur [0:1]  */
float inv_F_repartition_propagateur(float y, float mu, float sigma)
{
    return -mu * log(y) / sigma;
}

/* Variable aléatoire suivant la densité de probabilité
   du libre parcours d'un neutron */
part_trajectory_t one_traj(float mu, float sigma, source_t stype)
{
    float src = source(stype);
    return part_trajectory_t( src, inv_F_repartition_propagateur(ranf(), mu, sigma) + src);
}

// Echantillon de taille n de cette variable aléatoire
vector<part_trajectory_t> trajs(float mu, float sigma, int n, source_t stype)
{
    vector<part_trajectory_t> parts(n);
    for(int i = 0 ; i < n ; i++)
    {
	parts.at(i) = one_traj(mu, sigma, stype);
    }
    return parts;
}

/* Flux de particules dénombrant cet échantillon
   sur chaque segment de l'intervalle total */
distrib_t denom(const vector<part_trajectory_t>& parts, int nb_segs, float sigma)
{
    distrib_t distrib(nb_segs);
    for (auto const & val : parts) // sans indice
    {
	if (val.end>=1) distrib.above++; 
	else if (val.end< 0) distrib.below++;
	else distrib.segs.at( floor(nb_segs*val.end) )++;
    }
    distrib.above = distrib.above * nb_segs / (sigma*parts.size());
    distrib.segs  = distrib.segs  * nb_segs / (sigma*parts.size());
    distrib.below = distrib.below * nb_segs / (sigma*parts.size());
    return distrib;
}


int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc<7) {
	cout << "Enter mu, sigma, nb_particules, nb_segments, "
	     << "type source (constante=1, delta(0) = 2)"
	     << "output style (none=1, plot=2, file=3, all=4) [, filename]"
	     << endl;
	return EXIT_FAILURE;
    }
    float mu = atof(argv[1]),              // mu du problème
	sigma = atof(argv[2]);             // sigma du problème
    int nb_parts = floor(atof(argv[3])),   // Nombre de particules de l'échantillon (floor(atof()) pour pouvoir utiliser 1e6 etc...
	nb_segs  = floor(atof(argv[4])),   // Finesse de segmentation de l'intervalle pour calculer le flux
	stype_int = atoi(argv[5]),         // Source uniforme ou seulement entrante à gauche
	output_style = atoi(argv[6]);
    if ( !mu ||
	 !sigma ||
	 !nb_parts ||
	 !nb_segs ||
	 !stype_int ||
	 !output_style ){
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }
    source_t stype;
    if (stype_int==1) stype = CSTE;
    else stype = DELTA;

    string fname;
    if (output_style>2) {
	if (argc>7)
	    fname = argv[7];
	else
	    fname = "output_TP1";
    }
    
    // Courbe théorique
    vectorV2<float> X = linspace(0, 1, nb_segs), // vecteurs des abscisses
	Py (nb_segs);                            // vecteurs de la probabilité théorique
    for (int i=0; i<nb_segs; i++)
	Py.at(i) = Phi(X.at(i), mu, sigma, stype);

    // Monte Carlo
    vector<part_trajectory_t> parts = trajs(mu, sigma, nb_parts, stype);  // calcule des trajectoires 
    distrib_t distrib = denom(parts, nb_segs, sigma);                  // repartition dans les segments 

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

