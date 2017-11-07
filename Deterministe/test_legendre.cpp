#include "utils.hpp"
#include "vector_v2.hpp"
#include "legendre_polynomial.hpp"

using namespace std;
typedef vectorV2<double> v2d;

double f(double x) {
    return pow(x, 20);
}

v2d f_v2d(const v2d& X) {
    v2d Y (X.size());
    for (int i = 0; i < X.size(); ++i)
	Y.at(i) = f(X.at(i));
    return Y;
}

int main(int argc, char **argv) {
    int size = 100;
    v2d t (size),
	wts (size);
    p_quadrature_rule(size, t.data(), wts.data());
    v2d Y = f_v2d(t);
    cout << integ_gauss_legendre(Y, wts) << endl;
}
