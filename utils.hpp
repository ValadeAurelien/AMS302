#include <vector>http://timmurphy.org/2014/04/26/using-fork-in-cc-a-minimum-working-example/
#include <stdlib.h>
#include <string>
#include <iostream>
#include <cmath>
#include "vector_v2.hpp"

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
std::vector<float> linspace(float min, float max, unsigned nb_pts)
{
    std::vector<float> X(nb_pts);
    for (int i=0; i<nb_pts; i++) X[i] = min + i * (max-min) / nb_pts;
    return X;
}

/* 
** Fonctions de plot de 1 ou plusieurs vectors
*/
template<typename Tx, typename Ty>
int plot(const std::vector<Tx>& X, const std::vector<Ty>& Y, const std::string& opt = std::string(), const std::string& title = std::string())
{
    if (X.size() != Y.size()) {
	std::cout << "sizes don't match !" << std::endl;
	return -1;
    }
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ("set title '" + title + "' ; plot '-' " + opt + "\n").c_str());
    for (int i=0; i<X.size(); i++)
	fprintf(pipe, "%f %f\n", (float) X.at(i), (float) Y.at(i));
    fprintf(pipe, "pause -1;");
    fflush(pipe);
    return pclose(pipe);
}

template<typename Ty>
int plot(const std::vector<Ty>& Y, const std::string& opt = std::string(), const std::string& title = std::string())
{
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ("set title '" + title + "' ; plot '-' " + opt + "\n").c_str());
    for (int i=0; i<Y.size(); i++)
	fprintf(pipe, "%f\n", (float) Y.at(i));
    fflush(pipe);
    return pclose(pipe);
}
