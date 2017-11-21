#include "fe1D_solver.hpp"
#include <iostream>
#include <cmath>
#include <ctime>
#include <stdlib.h>
using namespace std;

int main(int argc, char **argv) {
    int size=atoi(argv[1]);
    double Mc = atof(argv[2]),
	Kc = atof(argv[3]);
    FE1DSolver_t solver(size, Mc, Kc);
    VectorXd L (size),
	X = linspacePts(0, 1, size), 
	Y(size);
    double dx = 1./(size-1);
    for (int i=0; i<size; i++) L(i) = 10*pow(i*dx-1/2.,2); //sin(M_PI*i*dx);
    solver.init();
    solver.setRightHand(L);
    clock_t c = clock();
    solver.solve(Y);
    clock_t d = clock();
    plot(X, Y, "w l t 'Y'", "");
    plot(X, L, "w l t 'L'", "");
    for (int i=0; i<size; ++i)
	cout << setw(15) << X(i) << setw(15) << L(i) << setw(15) << Y(i) << endl;
    cout << ((double) (d-c))/CLOCKS_PER_SEC << " " << (L-Y).norm()/L.norm() << endl;
    return 0;
}
