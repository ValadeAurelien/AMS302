#include <vector>
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
    part_trajectory_t (float _b, float _e, unsigned _n) : begin(_b), end(_e), nb_jumps(_n) {};

    float begin,
	end;
    unsigned nb_jumps;
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

/* Fonction qui renvoie la distance en fonction de la 
   proba uniforme sur [0:1]  */
float inv_F_repartition_propagateur(float y, float mu, float sigma_t)
{
    return -mu * log(y) / sigma_t;
}

/* Variable aléatoire suivant la densité de probabilité
   du libre parcours d'un neutron */
part_trajectory_t one_traj(float mu, float sigma_a, float sigma_s, source_t stype)
{
    float src = source(stype), 
	p_abs = ranf(),
	p_diff,
	sigma_t = sigma_a + sigma_s,
	sast = sigma_a/sigma_t,
	ssst = sigma_s/sigma_t,
	x = src;
    int n_max = 1000000,
	n=0;
    while ( p_abs>sast ) {
	n++;
	if ( n>n_max ) throw runtime_error("Particule did not finish");
	p_abs = ranf();
	p_diff = ranf();
	if ( p_diff>ssst )
	    mu = ranf();
	x += inv_F_repartition_propagateur(ranf(), mu, sigma_t);
    }
    return part_trajectory_t( src, x, n );
}

// Echantillon de taille n de cette variable aléatoire
vector<part_trajectory_t> trajs(float mu, float sigma_a, float sigma_s, int n, source_t stype)
{
    vector<part_trajectory_t> parts(n);
    for(int i = 0 ; i < n ; i++)
    {
	parts.at(i) = one_traj(mu, sigma_a, sigma_s, stype);
    }
    return parts;
}

/* Flux de particules dénombrant cet échantillon
   sur chaque segment de l'intervalle total */
distrib_t denom(const vector<part_trajectory_t>& parts, int nb_segs, float mu)
{
    distrib_t distrib(nb_segs);
    int seg_beg,
	seg_end;
    for (auto const & val : parts) {
	if (val.end>=1) {
	    distrib.above++;
	    seg_beg = ceil(nb_segs*val.begin);
	    for (int i=seg_beg; i<nb_segs; i++) distrib.segs.at(i)++;
	}
	else if (val.end<0) {
	    distrib.below++;
	    seg_beg = ceil(nb_segs*val.begin);
	    for (int i=0; i<seg_beg; i++) distrib.segs.at(i)++;
	}
	else {
	    seg_beg = ceil(nb_segs*val.begin);
	    seg_end = ceil(nb_segs*val.end);
	    if (seg_end > seg_beg)
		for (int i=seg_beg; i<seg_end; i++) distrib.segs.at(i)++;
	    else
		for (int i=seg_end; i<seg_beg; i++) distrib.segs.at(i)++;
	}
    }
    distrib.above = distrib.above / (abs(mu)*parts.size());
    distrib.segs  = distrib.segs  / (abs(mu)*parts.size());
    distrib.below = distrib.below / (abs(mu)*parts.size());    
    return distrib;
}

distrib_t denom_nb_jumps(const vector<part_trajectory_t>& parts, int max_nb_jumps)
{
    distrib_t distrib(max_nb_jumps+1);
    for (auto const& val : parts){
	if (val.nb_jumps > max_nb_jumps)
	    distrib.above++;
	else distrib.segs.at(val.nb_jumps)++;
    }
    return distrib;
}

int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc!=8) {
	cout << "Enter mu, sigma_a, sigma_s, nb_particules, nb_segments, type source (constante=1, delta(0) = 2), max_nb_jumps" << endl;
	return EXIT_FAILURE;
    }
    float mu = atof(argv[1]),              // mu du problème
	sigma_a = atof(argv[2]),
	sigma_s = atof(argv[3]);           // sigmas du problème
    int nb_parts = floor(atof(argv[4])),   // Nombre de particules de l'échantillon (floor(atof()) pour pouvoir utiliser 1e6 etc...
	nb_segs  = floor(atof(argv[5])),   // Finesse de segmentation de l'intervalle pour calculer le flux
	stype_int = atoi(argv[6]),         // Source uniforme ou seulement entrante à gauche
	max_nb_jumps = atoi(argv[7]);
    if ( !mu ||
	 !sigma_a ||
	 !sigma_s ||
	 !nb_parts ||
	 !nb_segs ||
	 !stype_int ||
	 !max_nb_jumps ){
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }
    source_t stype;
    if (stype_int==1) stype = CSTE;
    else stype = DELTA;

    vectorV2<float> X = linspace(0, 1, nb_segs);

    // Monte Carlo
    vector<part_trajectory_t> parts = trajs(mu, sigma_a, sigma_s, nb_parts, stype);  // calcule des trajectoires 
    distrib_t distrib = denom(parts, nb_segs, mu),                  // repartition dans les segments 
	distrib_nb_jumps = denom_nb_jumps(parts, max_nb_jumps);
    
    // Affichage
    for (int i=0; i<nb_segs; i++)
    	cout << X.at(i) << " " << distrib.segs.at(i) << endl;
    for (int i=0; i<=max_nb_jumps; i++)
    	cout << distrib_nb_jumps.segs.at(i) << endl;
    plot(X, distrib.segs, "w l", "set title 'MC'; set yrange [0:]; set ylabel '{/Symbol F}'; ");
    plot(distrib_nb_jumps.segs, "", "set title 'Nombre de sauts par particule'; set ylabel 'n';");
    cout << "#(above 1) / (below 0) : " << distrib.above << " / " << distrib.below << endl;
    return EXIT_SUCCESS;
}

