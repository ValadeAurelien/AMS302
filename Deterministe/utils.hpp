#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <cmath>
#include "vector_v2.hpp"
#include "legendre_polynomial.hpp"

/* 
** La structure distrib regroupe les infos nécéssaire 
** à la structuration d'une distributions des particules
** en dedans et en dehors du segment [0:1]
*/
struct distrib_t
{
    distrib_t () {}
    distrib_t (unsigned size) : segs(size, 0) {}
    distrib_t (const distrib_t& d) : below(d.below), above(d.above), segs(d.segs) {}
    float below,
	above;
    vectorV2<float> segs;
};

/* 
** Fonction qui renvoie un entier entre 0 et 1
*/
float ranf()
{
  return (float) rand()/RAND_MAX; //fdsqfaze
}

/*
** Analogue de la fonction de python (segmentation en n pts 
** du segment réel entre min et max)
*/
template <typename T>
vectorV2<T> linspace(float min, float max, unsigned nb_pts)
{
    std::vector<T> X(nb_pts);
    for (int i=0; i<nb_pts; i++) X[i] = (T) (min + i * (max-min) / nb_pts);
    return X;
}

/* 
** Fonctions de plot de 1 ou plusieurs vectors
*/
template<typename Tx, typename Ty>
int plot(const std::vector<Tx>& X, const std::vector<Ty>& Y, const std::string& plot_opt = std::string(), const std::string& pre_cmd = std::string())
{
    if (X.size() != Y.size()) {
	std::cout << "sizes don't match !" << std::endl;
	return -1;
    }
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ( pre_cmd + "plot '-' " + plot_opt + "\n" ).c_str());
    for (int i=0; i<X.size(); i++)
	fprintf(pipe, "%f %f\n", (float) X.at(i), (float) Y.at(i));
    fprintf(pipe, "pause -1;");
    fflush(pipe);
    return pclose(pipe);
}

template<typename Ty>
int plot(const std::vector<Ty>& Y, const std::string& plot_opt = std::string(), const std::string& pre_cmd = std::string())
{
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ( pre_cmd + "plot '-' " + plot_opt + "\n" ).c_str());
    for (int i=0; i<Y.size(); i++)
	fprintf(pipe, "%f\n", (float) Y.at(i));
    fflush(pipe);
    return pclose(pipe);
}


/*
** Fonction d'intégration d'une fonction sur un segment
*/
template<typename T>
double integ_gauss_legendre(const vectorV2<T>& Y,
			    const vectorV2<double>& Weights) {
    if ( Y.size() != Weights.size() )
	throw size_exception();

    double sum = 0;
    int size = Y.size();
    for (int i = 0; i < size; ++i)
	sum += Y.at(i)*Weights.at(i);
    return sum;
}
#endif
