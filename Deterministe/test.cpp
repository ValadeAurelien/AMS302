#include "fe1D_solver.hpp"
#include <iostream>
#include <cmath>
#include <ctime>
#include <stdlib.h>
using namespace std;

int main(int argc, char **argv) {
    int size=atoi(argv[1]);
    FE1DSolver_t solver(size, -1);
    VectorXd L (size);
    for (int i=0; i<size; i++) L(i) = sin((2*M_PI*i)/size);
    solver.setRightHand(L);
    VectorXd X;
    clock_t c = clock();
    solver.solve(X);
    clock_t d = clock();
    cout << ((double) (d-c))/CLOCKS_PER_SEC << " " << (L-X).norm()/size << endl;
    return 0;
}
