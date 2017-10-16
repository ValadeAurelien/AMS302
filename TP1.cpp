#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>
using namespace std;

// Variable aléatoire uniforme entre 0 et 1
float ranf()
{
  return (float) rand()/RAND_MAX;
}

/* Terme source, suivant qu'il est uniforme (typ_source vrai)
                  ou non-nul seulement en 0 (typ_source faux) */
float source(bool typ_source)
{
  if (typ_source)
    return ranf();
  else
    return 0;
}

/* Variable aléatoire suivant la densité de probabilité

du libre parcours d'un neutron */
float f_gene(float mu, float sigma, bool typ_source)
{ 

  return -mu * log(ranf()) / sigma + source(typ_source);
}

// Echantillon de taille n de cette variable aléatoire
vector<float> trajs(float mu, float sigma, int n, bool typ_source)
{
  vector<float> part(n);
  for(int i = 0 ; i < n ; i++)
  {
    part.at(i) = f_gene(mu,sigma);
  }
  return part;
}

/* Flux de particules dénombrant cet échantillon
   sur chaque segment de l'intervalle total */
vector<int> denom(const vector<float>& part, int nb_seg, bool typ_source)
{
  vector<int> flux(nb_seg + 2,0);
  for (auto const & val : part)
  {
    if (val>=1) flux.at(nb_seg+1)++; // Les deux cas où le particule ...
    else if (val< 0) flux.at(0)++;   // ... sort de l'intervalle
    else flux.at( floor(nb_seg*val) )++;
  }
  return flux;
}

int main(int argc, char**argv)
{
  if (argc!=5) return EXIT_FAILURE;
  float mu,
        sigma;
  int nb_part,
      nb_seg;
  bool typ_source;
  mu = atof(argv[1]);
  sigma = atof(argv[2]);
  nb_part = atoi(argv[3]);    // Nombre de particules de l'échantillon
  nb_seg  = atoi(argv[4]);    // Finesse de segmentation de l'intervalle pour calculer le flux
  typ_source = atoi(argv[5]); // Source uniforme ou seulement entrante à gauche
  
  vector<float> part = trajs(mu, sigma, nb_part,typ_source);
  //for(auto const & val : part)
  //  cout << val << endl;
  vector<int> flux = denom(part, nb_seg,typ_source);
  // for(int i = 1 ; i < nb_seg ; i++)
  //   cout << (float) i / nb_seg << " " << (float) flux.at(i) / nb_part << endl;
  return EXIT_SUCCESS;
}
