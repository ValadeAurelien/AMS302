#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <cmath>
#include <Eigen/Eigen>
#include "legendre_polynomial.hpp"

using namespace std;
using namespace Eigen;

//typedef Matrix<double, Dynamic, 1> VectorXd;
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
    VectorXd segs;
};

struct size_exception : std::exception
{
    const char* what() const noexcept { return "Size(s) do(es) not match !"; }
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
VectorXd linspaceSegs(float min, float max, Index nb_segs)
{
    VectorXd X (nb_segs);
    double dx = (max-min) / nb_segs;
     for (int i=0; i<nb_segs; i++)
     	X(i) = min + i * dx ;
     return X;
}

/*
** Analogue de la fonction de python (segmentation en n pts 
** du segment réel entre min et max)
*/
VectorXd linspacePts(float min, float max, Index nb_pts)
{
    VectorXd X (nb_pts);
    double dx = (max-min) / (nb_pts-1);
     for (int i=0; i<nb_pts; i++)
     	X(i) = min + i * dx ;
     return X;
}
/* 
** Fonctions de plot de 1 ou plusieurs vectors
*/
int plot(const VectorXd& X, const VectorXd& Y, const string& plot_opt = string(), const string& pre_cmd = string())
{
    if (X.size() != Y.size()) {
	cout << "sizes don't match !" << endl;
	return -1;
    }
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ( pre_cmd + "plot '-' " + plot_opt + "\n" ).c_str());
    for (int i=0; i<X.size(); i++)
	fprintf(pipe, "%f %f\n", (float) X(i), (float) Y(i));
    fprintf(pipe, "pause -1;");
    fflush(pipe);
    return pclose(pipe);
}

int plot(const VectorXd& Y, const string& plot_opt = string(), const string& pre_cmd = string())
{
    FILE *pipe = popen("gnuplot -p", "w");
    fprintf(pipe, ( pre_cmd + "plot '-' " + plot_opt + "\n" ).c_str());
    for (int i=0; i<Y.size(); i++)
	fprintf(pipe, "%f\n", (float) Y(i));
    fflush(pipe);
    return pclose(pipe);
}

/*
** Fonctions d'intégration d'une fonction sur un segment
*/
void get_gauss_lengendre_pts_wghts(VectorXd& X, VectorXd& W) {
    p_quadrature_rule(X.rows(), X.data(), W.data());
}

void get_cst_pts_wghts(VectorXd& X, VectorXd&W) {
    int size = X.rows();
    double dx = 1./(size-1);
    X = linspacePts(0, 1, size);
    for (int i=0; i<size; ++i) W(i) = 2*dx;
}

double integ_gauss_legendre(const VectorXd& Y,
			    const VectorXd& Weights) {
    if ( Y.rows() != Weights.rows() )
	throw size_exception();

    double sum = 0;
    int size = Y.size();
    for (int i = 0; i < size; ++i)
	sum += Y(i)*Weights(i);
    return sum;
}
#endif
