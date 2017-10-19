#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
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
typedef vector<float> part_trajectory_t;

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
void one_traj(part_trajectory_t& part, float sigma_a, float sigma_s, source_t stype)
{
    float src = source(stype), 
	p_abs = ranf(),
	p_diff,
	sigma_t = sigma_a + sigma_s,
	sast = sigma_a/sigma_t,
	x = src,
	mu = 2*ranf()-1;
    int n_max = 1000000,
	n=0;
    while ( p_abs>sast ) {
	if ( n++>n_max )
	    throw runtime_error("Particule did not finish");
	if ( x<0 || x>1 ) break;
	p_abs = ranf();
	p_diff = ranf();
	if ( p_diff<1-sast ) {
	    part.push_back(x);
	    mu = 2*ranf()-1;
	}
	x += inv_F_repartition_propagateur(ranf(), mu, sigma_t);
    }
}

// Echantillon de taille n de cette variable aléatoire
void trajs(vector<part_trajectory_t>& parts, float sigma_a, float sigma_s, source_t stype)
{
    for(auto& p : parts) 
	one_traj(p, sigma_a, sigma_s, stype);
}

/* Flux de particules dénombrant cet échantillon
   sur chaque segment de l'intervalle total */   

void denom(distrib_t& distrib, const vector<part_trajectory_t>& parts, int nb_segs, float sigma_t)
{
    for (auto const & part : parts) // sans indice
	for (auto const& pt : part) {
	    if (pt>=1) distrib.above++; 
	    else if (pt<0) distrib.below++;
	    else distrib.segs.at( floor(nb_segs*pt) )++;
	}
    distrib.above = distrib.above * nb_segs / (sigma_t*parts.size()); // on divise par l'intégrale du propagateur
    distrib.segs  = distrib.segs  * nb_segs / (sigma_t*parts.size());
    distrib.below = distrib.below * nb_segs / (sigma_t*parts.size());
}

void denom_nb_jumps(distrib_t& distrib, const vector<part_trajectory_t>& parts, int max_nb_jumps)
{
    for (auto const& part : parts){
	if (part.size() > max_nb_jumps)
	    distrib.above++;
	else distrib.segs.at(part.size())++;
    }
}

int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc<8) {
	cout << "Enter sigma_a, sigma_s, nb_particules, nb_segments, "
	     << "type source (constante=1, delta(0) = 2), max_nb_jumps, "
	     << "output style (none=1, plot=2, file=3, all=4) [, filename]"
	     << endl;
	return EXIT_FAILURE;
    }
    float sigma_a = atof(argv[1]),
	  sigma_s = atof(argv[2]);           // sigmas du problème
    int nb_parts = floor(atof(argv[3])),   // Nombre de particules de l'échantillon (floor(atof()) pour pouvoir utiliser 1e6 etc...
	nb_segs  = floor(atof(argv[4])),   // Finesse de segmentation de l'intervalle pour calculer le flux
	stype_int = atoi(argv[5]),         // Source uniforme ou seulement entrante à gauche
	max_nb_jumps = atoi(argv[6]),
	output_style = atoi(argv[7]);
    if ( !sigma_s ||
	 !nb_parts ||
	 !nb_segs ||
	 !stype_int ||
	 !max_nb_jumps ||
	 !output_style ){
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }
    source_t stype;
    if (stype_int==1) stype = CSTE;
    else stype = DELTA;
    string fname;
    if (output_style>2) {
	if (argc>8)
	    fname = argv[8];
	else
	    fname = "output_TP1_diff";
    }
    
    // X pour plot
    vectorV2<float> X = linspace(0, 1, nb_segs);

    // Monte Carlo
    vector<part_trajectory_t> parts (nb_parts);
    distrib_t distrib (nb_segs),
	distrib_nb_jumps (max_nb_jumps+1);
    trajs(parts, sigma_a, sigma_s, stype);  // calcule des trajectoires 
    denom(distrib, parts, nb_segs, sigma_a+sigma_s);             // repartition dans les segments 
    denom_nb_jumps(distrib_nb_jumps, parts, max_nb_jumps);
    
    // Affichage
    if (output_style>2) {
	ofstream file (fname, fstream::out);
	file << "#";
	for (int i=0; i<argc; i++) file << argv[i] << " ";
	file << endl;
	file << "#distrib" << endl << "#X     Phi" << endl;
	for (int i=0; i<nb_segs; i++)
	    file << X.at(i) << " " << distrib.segs.at(i) << endl;
	file << "#(above 1) / (below 0) : "
	     << distrib.above << " / " << distrib.below << endl;
	file << "#nb_jumps" << endl << "#nj    npart" << endl; 
	for (int i=0; i<=max_nb_jumps; i++)
	    file << i << " " << distrib_nb_jumps.segs.at(i) << endl;
	file.close();
    }
    if (!(output_style % 2)) {
	plot(X, distrib.segs, "w l", "set title 'MC'; set yrange [0:]; set ylabel '{/Symbol F}'; ");
	plot(distrib_nb_jumps.segs, "", "set title 'Nombre de sauts par particule'; set ylabel 'n';");
	cout << "#(above 1) / (below 0) : " << distrib.above << " / " << distrib.below << endl;
    }
    return EXIT_SUCCESS;
}

