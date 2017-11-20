#ifndef FE1D_PROBLEM_HEADER
#define FE1D_PROBLEM_HEADER
#include <Eigen/Eigen>
#include <iostream>
#include "utils.hpp"
using namespace Eigen;
using namespace std;

struct FE1DSolver_t {
    double K_coef;
    SparseMatrix<double> M, K, A;
    ConjugateGradient<SparseMatrix<double> > solver;
    VectorXd B;
    int size;

    FE1DSolver_t(int _size, double _Kc)
	: size(_size), M((Index) _size, (Index) _size), K((Index) _size, (Index) _size),
	  K_coef (_Kc), B(_size) {
	double dx = (double) 1/(size-1);
	M.insert(0, 0) = dx   ; M.insert(0, 1) = dx/2 ;
	K.insert(0, 0) = 2*dx ; K.insert(0, 1) = -dx  ;
	for (int i=1; i<size-1; ++i) {
	    M.insert(i, i-1) = dx/2 ; M.insert(i, i) = dx   ; M.insert(i, i+1) = dx/2 ;
	    K.insert(i, i-1) = -dx  ; K.insert(i, i) = 2*dx ; K.insert(i, i+1) = -dx  ;
	}
	A = M + K_coef * K;
	solver.compute(A);
	if(solver.info()!=Success) {
	    cout << "Impossible to compute this matrix"<< endl;
	    return;
	}
	
    }
    void setRightHand(VectorXd& L) {
	B = M * L;
    }
    void solve(VectorXd& X) {
	X = solver.solve(B);
	if(solver.info()!=Success) {
	    cout << "Impossible to solve the system"<< endl;
	    //throw exception("Impossible to solve the system");
	    return;
	}
    }
};

#endif
