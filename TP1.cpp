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
    case DELTA:
	return 1/mu*exp(-sigma/mu*x);
    }
}

/* Fonction qui renvoie la distance en fonction de la 
   proba uniforme sur [0:1]  */
float inv_F_repartition_propagateur(float x, float mu, float sigma)
{
    return -mu * log(x) / sigma;
}

/* Variable aléatoire suivant la densité de probabilité
   du libre parcours d'un neutron */
float f_gene(float mu, float sigma, source_t stype)
{ 
    return inv_F_repartition_propagateur(ranf(), mu, sigma) + source(stype);
}

// Echantillon de taille n de cette variable aléatoire
vectorV2<float> trajs(float mu, float sigma, int n, source_t stype)
{
    vectorV2<float> part(n);
    for(int i = 0 ; i < n ; i++)
    {
	part.at(i) = f_gene(mu, sigma, stype);
    }
    return part;
}

/* Flux de particules dénombrant cet échantillon
   sur chaque segment de l'intervalle total */
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
    distrib.segs  = distrib.segs  * nb_segs / parts.size();
    distrib.below = distrib.below * nb_segs / parts.size();
    return distrib;
}

int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc!=6) {
	cout << "Enter mu, sigma, nb_particules, nb_segments, type source (constante=1, delta(0) = 2)" << endl;
	return EXIT_FAILURE;
    }
    float mu = atof(argv[1]),              // mu du problème
	sigma = atof(argv[2]);             // sigma du problème
    int nb_parts = floor(atof(argv[3])),   // Nombre de particules de l'échantillon (floor(atof()) pour pouvoir utiliser 1e6 etc...
	nb_segs  = floor(atof(argv[4])),   // Finesse de segmentation de l'intervalle pour calculer le flux
	stype_int = atoi(argv[5]);         // Source uniforme ou seulement entrante à gauche
    if ( !mu ||
	 !sigma ||
	 !nb_parts ||
	 !nb_segs ||
	 !stype_int ){
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }
    source_t stype;
    if (stype_int==1) stype = CSTE;
    else stype = DELTA;

    vectorV2<float> X = linspace(0, 1, nb_segs), // vecteurs des abscisses
	Py (nb_segs);                            // vecteurs de la probabilité théorique
    for (int i=0; i<nb_segs; i++)
	Py.at(i) = Phi(X.at(i), mu, sigma, stype);
    
    vectorV2<float> parts = trajs(mu, sigma, nb_parts, stype);  // calcule des trajectoires 
    distrib_t distrib = denom(parts, nb_segs);                  // repartition dans les segments 
    plot(X, distrib.segs, "w l");
    plot(X, Py, "w l");
    cout << "#(above 1) / (below 0) : " << distrib.above << " / " << distrib.below << endl
	 << "#diff normalized : " << (distrib.segs-Py).norm()/Py.norm() << endl;
    return EXIT_SUCCESS;
}
