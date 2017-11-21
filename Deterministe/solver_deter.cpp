#include <Eigen/Eigen>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <stdexcept>
#include "utils.hpp"
#include "fe1D_solver.hpp"

#define WIDTH 15
#define PRECISION 5
#define MAX_COUNT 10000
//using namespace std;
//using namespace Eigen;

// Deux types de source
enum func_t {
    CSTE = 1,
    DELTA = 2,
    TWO_STEPS = 3
};

// Masques pour la sortie
enum output_t {
    OUTPUT_PLOT = 1,
    OUTPUT_FILE = 2
};

// Plusieurs types d'experience
enum expe_type_t {
    NO_SCATTERING = 1,
    SCATTERING = 2,
    SCATTERING_DSA = 3
};

/* Terme source */
struct source_func_t {
    int sourcet;
    double arg1;

    source_func_t(int _st, double _arg1)
	: sourcet (_st), arg1(_arg1) {}

    double operator()(const double& x) const {
	switch(sourcet) {
	case CSTE :
	    return arg1;
	case DELTA:
	    return (x==0)*arg1;
	default:
	    return 0;
	}
    }
};

/* Coefficient de l'équation différentielle */
struct sigma_func_t {
    int sigmat;
    double arg1;

    sigma_func_t(int _sigmat, double _arg1)
	: sigmat(_sigmat), arg1(_arg1) {}

    double operator()(const double& x) const {
	switch(sigmat) {
	case CSTE:
	    return arg1;
	case TWO_STEPS:
	    if (x<0.3)
		return arg1;
	    if (x<0.7)
		return 3*arg1;
	    else
		return arg1;
	    break;
	default:
	    return 0;
	}
    }
};

/* Calcul de la solution exacte 
   beaucoup de cas particuliers */
struct phi_func_t {
    double mu;
    sigma_func_t& sigma_a_func,
	&sigma_s_func;
    source_func_t& source_func;
    
    phi_func_t(double _mu, sigma_func_t& _siaf, sigma_func_t& _sisf, source_func_t& _sof)
	: mu(_mu), sigma_a_func(_siaf), sigma_s_func(_sisf), source_func(_sof) {}

    double sigma_t(const double& x) const {
	return sigma_a_func(x)+sigma_s_func(x);
    }
    
    double operator()(const double& x) const {
	switch (sigma_a_func.sigmat) {
	case CSTE: // type de sigma
	    switch(source_func.sourcet) {
	    case CSTE :  // type de source
		if (mu>0)
		    return 1/sigma_t(x)*(1-exp(-sigma_t(x)/mu*x));
		else
		    return 1/sigma_t(x)*(1-exp(-sigma_t(x)/mu*(x-1)));
	    case DELTA :
		if (mu>0)
		    return 1/mu*exp(-sigma_t(x)/mu*x);
		else
		    return 0;
	    default:
		return 0;
	    }
	    break;
	case TWO_STEPS: // type de sigma
	    switch(source_func.sourcet) {
	    case DELTA : // type de source
		if (mu>0)
		    if (x<0.3) return 1/mu*exp(-x/mu);
		    else if (x<0.7) return 1/mu*exp((0.6-3*x)/mu);
		    else return 1/mu*exp(-(x+0.8)/mu);
		else
		    return 0;
	    default:
		return 0;
	    }
	    break;
	}
    }
};

/* Wrapper du solver FE 1D du fichier fe1D_solver.hpp */
struct DSA_solver_t {
    const sigma_func_t& sigma_a_func,
	&sigma_s_func;
    VectorXd& dQ,
	&f, sdQ;
    FE1DSolver_t FE1DSolver;

    DSA_solver_t(const sigma_func_t& _saf, const sigma_func_t& _ssf,
		 VectorXd& _dQ, VectorXd& _f)
	: sigma_a_func(_saf), sigma_s_func(_ssf), dQ(_dQ), f(_f),
	  FE1DSolver ( f.rows(), sigma_a_func(0),
		       1/(3.*sigma_s_func(0)+sigma_a_func(0)) )	{
    }

    void init() {
	FE1DSolver.init();
    }

    int solve() {
	sdQ = sigma_s_func(0)*dQ;
	FE1DSolver.setRightHand(sdQ);
	FE1DSolver.solve(f);
    }
};		     

/* Phi déterministe avec diffusion 
   tout est la dessous... */
VectorXd Phi_deter_many_mu(Index nb_pts, int nb_pts_mu,
			   sigma_func_t& sigma_a_func, sigma_func_t& sigma_s_func,
			   source_func_t& source_func,
			   double epsilon, bool DSA) {
    // Initialisation
    VectorXd X = linspacePts(0,1,nb_pts),
	Q = X.unaryExpr(source_func),
	Qt = Q,
	dQ = Q-Qt,
	phi (nb_pts),
	phit (nb_pts),
	f_DSA (nb_pts),
	f_DSAt (nb_pts),
	vec_mus (nb_pts_mu),
	weights (nb_pts_mu);
    double phip, phim, eta_m, eta_p,
	x = 0,
	dx = (double) 1/(nb_pts-1);
    int count = 0,
	nmu;
    DSA_solver_t DSA_Solver (sigma_a_func, sigma_s_func, dQ, f_DSA);
    if (DSA)
	DSA_Solver.init();

    // Calcul des points d'interpolation des mu.
    // La fonction est définie dans utils.hpp
    get_gauss_lengendre_pts_wghts(vec_mus, weights);
    //get_cst_pts_wghts(vec_mus, weights);
    clock_t timer = clock();
    do {
	// Réinitialisation de phi
	phi.fill(0);
	// Si on a fait trop de tours on s'en va !
	if (count++ > MAX_COUNT) {
	    cout << endl << "norm & eps : " << dQ.norm()
		 << " " << epsilon << " " << endl;
	    cout.flush();				      
	    throw range_error("Could not converge...");
	}
	// Tous les 100 tours on donne des nouvelles
	if (! (count%100) ) {
	    cout << "\r#Nb loops : " << setw(5) << count
		 << " et |dQ| = " << setw(7) << setprecision(4) << dQ.norm();
	    cout.flush();
	}
	// Mu négatifs (les premiers de l'interpolation)
	for (nmu=0; nmu<nb_pts_mu/2; ++nmu) {
	    switch (source_func.sourcet) {
	    case DELTA:
		continue;
		break;
	    case CSTE:
		phip = 0;                     
		break;
	    }
	    x=1;
	    for (int nx=nb_pts-1; nx>=0; --nx) {
		eta_m = 2*abs(vec_mus(nmu)) -
		    dx*(sigma_a_func(x)+sigma_s_func(x));
		eta_p = 2*abs(vec_mus(nmu)) +
		    dx*(sigma_a_func(x)+sigma_s_func(x));
		phim = (2*dx*Q(nx) + eta_m * phip) / eta_p;
		phi(nx) += 1/2.*weights(nmu)*(phip+phim);
		phip = phim;
		x -= dx;
	    }
	}
	// Mu positifs
	for (nmu=nmu; nmu<nb_pts_mu; ++nmu) {
	    if (vec_mus(nmu)==0)
		continue;
	    switch (source_func.sourcet) {
	    case DELTA: 
		phim=1/vec_mus(nmu);
		break;
	    case CSTE:
		phim=0;
		break;
	    }
	    x=0;
	    for (int nx=0; nx<nb_pts; ++nx) {
		eta_m = 2*abs(vec_mus(nmu)) -
		    dx*(sigma_a_func(x)+sigma_s_func(x));
		eta_p = 2*abs(vec_mus(nmu)) +
		    dx*(sigma_a_func(x)+sigma_s_func(x));
		phip = (2*dx*Q(nx) + eta_m * phim) / eta_p;
		phi(nx) += 1/2.*weights(nmu)*(phip+phim);
		phim = phip;
		x += dx;
	    }
	}
	// Mise a jour du terme source
	Qt = Q;
	Q = X.unaryExpr(source_func) + sigma_s_func(0)/2. * phi;
	dQ = Q-Qt;
	// Si DSA, on finit le boulot
	if (DSA) {
	    f_DSAt = f_DSA;
	    // On lance une exception si on a pas pu résoudre le
	    // système (ca n'arrive pas si le problème est bien
	    // posé... 
	    try { DSA_Solver.solve(); }
	    catch (no_convergence& e) {
		cout << "After " << count << " iterations" << endl;
		throw e;
	    }
	    Q += 1./2*(f_DSAt + f_DSA);
	    dQ = Q-Qt;
	}
    } while (dQ.norm()>epsilon);
    // On affiche quelques informations intéressantes et on s'en va.
    timer = clock() - timer;
    cout << endl << "#Exited loop after "
	 << setw(5) << count << " iterations in "
	 << setw(7) << setprecision(4) << (double) timer/CLOCKS_PER_SEC << "sec with |dQ| = "
	 << setw(7) << setprecision(4) << dQ.norm() << endl;
    return phi;
}

/* Fonction pour la résolution sans diffusion 
   quelques cas particuliers pour les conditions aux limites
   ainsi que pour le signe de mu */
VectorXd Phi_deter(int nb_pts, double mu,
		   sigma_func_t& sigma_t_func,
		   source_func_t& source_func) {
    VectorXd phi(nb_pts);
    double eta_p, eta_m,
	Stot, tau, x,
	dx = (double) 1 / nb_pts;
    if (mu>0) {
	switch (source_func.sourcet) {
	case DELTA: 
	    phi(0)=1/mu;
	    break;
	case CSTE:
	    phi(0)=0;
	    break;
	}
	for (int i=1; i<nb_pts; i++) {
	    x = i * dx;
	    eta_p = sigma_t_func(x) / 2. + mu * nb_pts;
	    eta_m = sigma_t_func(x) / 2. - mu * nb_pts;
	    Stot = (source_func(x)+source_func(x-dx))/2;
	    phi(i) = 1/eta_p * ( Stot - eta_m * phi(i-1) ) ;
	}
    }
    else {
	switch (source_func.sourcet) {
	case DELTA : 
	    for (int i=1; i<nb_pts; i++) { phi(i) = 0; }
	    break;
	case CSTE :
	    phi(nb_pts-1)= 0;                    
	    for (int i=nb_pts-2; i>=0; i--) {
		x = i * dx;
		eta_p = sigma_t_func(x) / 2. + mu * nb_pts;
		eta_m = sigma_t_func(x) / 2. - mu * nb_pts;
		Stot = (source_func(x)+source_func(x+dx))/2;
		phi(i) = 1/eta_m * ( Stot - eta_p * phi(i+1) ) ;
	    }
	    break;
	default :
	    throw invalid_argument("Mauvais argument de type d'expérience ou mauvais signe de mu");
	}
    }
    return phi;
}

/* Fonction principale... */
int main(int argc, char**argv) {
    // Vérification des arguments en entrée 
    if (argc<12) {
	cout << "Enter experimentation type (no_scattering=1, scattering=2), \n"
	     << "nb_segments, \n"
	     << "mu, nb_pts_mu, \n"
	     << "type sigma_a (constante=1), sa_arg1, \n"
	     << "type sigma_s (constante=1), ss_arg1, \n"
	     << "type source (constante=1, delta(0) = 2), s_arg1 \n"
	     << "epsilon CV, \n"
	     << "output style (none=0, plot=1, file=2, all=3) [, filename]"
	     << endl;
	return EXIT_FAILURE;
    }
    int    expe_type    = atoi(argv[1]);
    int    nb_segs      = floor(atof(argv[2]));
    double mu           = atof(argv[3]);            
    int    nb_pts_mu    = floor(atof(argv[4])); 
    int    sigma_at     = atoi(argv[5]);
    double sigma_a_arg1 = atof(argv[6]);
    int    sigma_st     = atoi(argv[7]);
    double sigma_s_arg1 = atof(argv[8]);
    double epsilon      = atof(argv[9]);
    int    sourcet      = atoi(argv[10]);
    double source_arg1  = atof(argv[11]);
    int    output_style = atoi(argv[12]);
    if ( !mu ||
	 !nb_segs ||
	 !nb_pts_mu ||
	 !sigma_at ||
	 !sigma_st ||
	 !sourcet ){
	throw invalid_argument("Mauvais arguments, ou mauvais types...");
    }    
    string fname;
    if (output_style & OUTPUT_FILE) {
	if (argc>13)
	    fname = argv[13];
	else
	    fname = "output";
    }

    // En fait on donne le nombre de segments, il y a un point en plus
    // ce point est en fait celui en x=1
    Index nb_pts = nb_segs+1;
    switch(expe_type) {
    case NO_SCATTERING:
	switch(sourcet) {
	case DELTA:
	    source_arg1=1/mu;
	}
	sigma_s_arg1=0;
	break;
    }
    // Initialisation des structes. Ca evite de passer tous
    // arguments aux autres fonctions derrière
    source_func_t source_func(sourcet,  source_arg1);
    sigma_func_t sigma_a_func(sigma_at, sigma_a_arg1);
    sigma_func_t sigma_s_func(sigma_st, sigma_s_arg1);
    phi_func_t phi_func(mu, sigma_a_func, sigma_s_func, source_func);

    // Vérification que le calcul est faisable avec
    // les données en entrées.
    double eta_m;
    switch(expe_type) {
    case NO_SCATTERING:
	eta_m = sigma_a_func(0.5)/2-abs(mu)*nb_pts;
    	if (eta_m>0) {
    	    cout << "eta_m positif : " << eta_m << endl;
    	    return EXIT_FAILURE;
    	}
    	break;
    case SCATTERING:
    case SCATTERING_DSA:
	eta_m = sigma_a_func(0.5)/2-2*1./nb_pts_mu*nb_pts;
    	if (eta_m>0) { //pas exactement la bonne formule avec la quadrature de gauss 
    	    cout << "eta_m positif : " << eta_m << endl;
    	    return EXIT_FAILURE;
    	}
    	break;
    }
    
    VectorXd X = linspacePts(0, 1, nb_pts);
    VectorXd Py (nb_pts),
	Py_deter(nb_pts);
    bool DSA=false;
    // On fait vraiment le travail ici ->
    switch(expe_type) {
    case NO_SCATTERING:
	Py  = X.unaryExpr(phi_func), // vecteurs de la probabilité théorique
	Py_deter = Phi_deter(nb_pts, mu, sigma_a_func, source_func);
	break;
    case SCATTERING_DSA:
	DSA = true;
    case SCATTERING:
	Py_deter = Phi_deter_many_mu(nb_pts, nb_pts_mu, sigma_a_func,   
				     sigma_s_func, source_func,
				     epsilon, DSA);
	break;
    default:
	throw invalid_argument("Mauvais argument de type d'expérience");
    }
    
    // Ecriture dans un fichier et plot (la fonction est dans utils.hpp) 
    if (output_style & OUTPUT_FILE) {
    	ofstream file (fname, fstream::out);
    	file << "#";
    	for (int i=0; i<argc; i++) file << argv[i] << " ";
    	file << endl;
    	file << "#distrib" << endl
	     << setw(WIDTH) << "#X"
	     << setw(WIDTH) << "Phi(APPROX)"
	     << setw(WIDTH) << "Phi(th)"
	     << endl;
    	for (int i=0; i<nb_pts; i++)
    	    file << setw(WIDTH) << setprecision(PRECISION) << X(i)
		 << setw(WIDTH) << setprecision(PRECISION) << Py_deter(i)
		 << setw(WIDTH) << setprecision(PRECISION) << Py(i) << endl;
    	file << "#diff normalized : " << (Py_deter-Py).norm()/Py.norm()
    	     << endl;
    	file.close();
    }
    if (output_style & OUTPUT_PLOT) {
    	plot(X, Py_deter, "w l", "set title 'Approx'; set yrange [0:]; ");
    	plot(X, Py, "w l", "set title 'théorique'; set yrange [0:]; ");
    	cout << "#diff normalized : " << (Py_deter-Py).norm()/Py.norm() << endl;
    }
    return EXIT_SUCCESS;
}

