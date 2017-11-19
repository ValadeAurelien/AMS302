#include <Eigen/Eigen>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdexcept>

using namespace std;
typedef Evector<typename T> Eigen::Matrix<T, Eigen::Dynamic, 1>;

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
    case CSTE :  // FAUX
	if (mu>0)
	    return 1/sigma*(1-exp(-sigma/mu*x));
	else
	    return 1/sigma*(1-exp(-sigma/mu*(x-1)));
    case DELTA :
        if (mu>0)
	    return 1/mu*exp(-sigma/mu*x);
        else
            return 0;
    }
}

/* Phi déterministe */
Evector<float> Phi_deter(int nb_segs, float mu, float sigma, source_t stype)
{
    Evector<float> vec_deter(nb_segs+1,0);
    switch(stype){
    case DELTA :
        if (mu>0) {
            vec_deter(0) = 1./mu;
            float k   = sigma / (nb_segs * mu);
            float tau = (1. - k/2) / (1. + k/2);
	    for (int i=1; i<nb_segs+1; i++)
                vec_deter(i) = vec_deter(i-1) * tau;
        }
        return vec_deter;
    }
}


int main(int argc, char**argv)
{
    // Vérification des arguments en entrée 
    if (argc<6) {
	cout << "Enter mu, sigma, nb_segments, "
	     << "type source (constante=1, delta(0) = 2)"
	     << "output style (none=1, plot=2, file=3, all=4) [, filename]"
	     << endl;
	return EXIT_FAILURE;
    }
    float mu = atof(argv[1]),              // mu du problème
	sigma = atof(argv[2]),             // sigma du problème
	nb_segs  = floor(atof(argv[3])),   // Finesse de segmentation de l'intervalle pour calculer le flux
	stype_int = atoi(argv[4]),         // Source uniforme ou seulement entrante à gauche
	output_style = atoi(argv[5]);
    if ( !mu ||
	 !sigma ||
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
	if (argc>6)
	    fname = argv[6];
	else
	    fname = "output_TP1";
    }

    // Courbe théorique
    Evector<float> X = linspace(0, 1, nb_segs+1), // vecteurs des abscisses
	Py (nb_segs+1);                            // vecteurs de la probabilité théorique
    for (int i=0; i<nb_segs+1; i++)
	Py.at(i) = Phi(X.at(i), mu, sigma, stype);
    // Courbe déterministe
    Evector<float> Py_deter = Phi_deter(nb_segs, mu, sigma, stype);

    // Affichage
    if (output_style>2) {
	ofstream file (fname, fstream::out);
	file << "#";
	for (int i=0; i<argc; i++) file << argv[i] << " ";
	file << endl;
	file << "#distrib" << endl << "#X     Phi(MC)     Phi(th)" << endl;
	for (int i=0; i<nb_segs; i++)
	    file << X(i) << " " << Py_deter(i) << " " << Py(i) << endl;
	file << "#diff normalized : " << (Py_deter-Py).norm()/Py.norm()
	     << endl;
	file.close();
    }
    if (output_style == 2 || output_style == 4) {
	plot(X, Py_deter, "w l", "set title 'MC'; set yrange [0:]; ");
	plot(X, Py, "w l", "set title 'théorique'; set yrange [0:]; ");
	cout << "#diff normalized : " << (Py_deter-Py).norm()/Py.norm() << endl;
    }
    return EXIT_SUCCESS;
}

