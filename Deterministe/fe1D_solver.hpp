#ifndef FE1D_PROBLEM_HEADER
#define FE1D_PROBLEM_HEADER
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include "utils.hpp"

using namespace Eigen;
using namespace std;

void my_print(const MatrixXd& A) {
    for (int i=0; i<A.rows(); i++) {
	for (int j=0; j<A.cols(); j++)
	    cout << setw(12) << setprecision(7) << A(i, j);
	cout << endl;
    }
}

struct no_convergence : exception {
    const char* what() const noexcept {
	return "Could not converge...";
    };
};

struct FE1DSolver_t {
    int size;
    double M_coef, K_coef;
    SparseMatrix<double> M, K, A, Id_mask;
    SimplicialLDLT<SparseMatrix<double> > solver;
    VectorXd X, B;

    FE1DSolver_t() = default;
    FE1DSolver_t(int _size, double _Mc, double _Kc)
	: size(_size), M((Index) size, (Index) size), K((Index) size, (Index) size),
	  Id_mask((Index) size-2, (Index) size), M_coef(_Mc), K_coef(_Kc), B(size), X(size) {
    }

    void init() {
	int s = size-1;
	double dx = (double) 1/s;

	M.insert(0, 0) = dx  ; M.insert(0, 1) = dx/2 ;
	K.insert(0, 0) = 2*s ; K.insert(0, 1) = -s   ;
	for (int i=1; i<size-1; ++i) {
	    Id_mask.insert(i-1, i) = 1;
	    M.insert(i, i-1) = dx/2 ; M.insert(i, i) = dx   ; M.insert(i, i+1) = dx/2 ;
	    K.insert(i, i-1) = -s   ; K.insert(i, i) = 2*s  ; K.insert(i, i+1) = -s   ;
	}
	M.insert(size-1, size-2) = dx  ; M.insert(size-1, size-1) = dx/2 ;
	K.insert(size-1, size-2) = 2*s ; K.insert(size-1, size-1) = -s   ;	
	A = Id_mask*(M_coef * M + 2*K_coef * K)*Id_mask.transpose();
	solver.compute(A);
	if(solver.info()!=Success) {
	    cout << "Impossible to compute this matrix"<< endl;
	    return;
	}
    }
    
    void setRightHand(VectorXd& L) {
	B = Id_mask * M * L;
    }
    void solve(VectorXd& Y) {
	X = solver.solve(B);
	Y = Id_mask.transpose()*X;
	if(solver.info()!=Success) {
	    cout << "Impossible to solve the system... " << endl
		 << "B (" << B.rows() << ") = " << endl << B << endl
		 << "X (" << X.rows() << ") = " << endl << X << endl;
	    throw no_convergence();
	}
	
    }
};

#endif
