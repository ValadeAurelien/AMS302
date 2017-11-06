#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <curand.h>
#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#define MIN(a, b) (((a)<(b)?(a):(b)))


__device__
float source(float rand_b) { return rand_b; }

__device__
float f_gene(float mu, float sigma, float rand_a, float rand_b)
{
    return -mu * logf( rand_a ) / sigma + source(rand_b);
}

__global__
void trajs(float mu, float sigma, float* parts, unsigned nb_parts, float* rands_a, float* rands_b )
{
    int x = threadIdx.x+blockIdx.x*blockDim.x;

    if (x>=nb_parts) return;
    
    parts[x] = f_gene(mu, sigma, rands_a[x], rands_b[x]);
}

__global__
void make_distrib(float* parts,
		  unsigned nb_parts,
		  unsigned* distrib,
		  unsigned nb_segs,
		  unsigned* below,
		  unsigned* above,
		  float min,
		  float max,
		  unsigned nb_threads)
{
    unsigned x = threadIdx.x+blockIdx.x*blockDim.x;

    if (x>=nb_threads) return;

    unsigned range_size = floorf((float) nb_parts/nb_threads),
	i = x*range_size;
    int seg = floorf( (float) (parts[i]-min)/(max-min)*nb_segs );
    for (i++; i<(x+1)*range_size; i++){
	if ( floorf( (float) (parts[i]-min)/(max-min)*nb_segs ) > seg )
	    seg = (int) floorf( (float) (parts[i]-min)/(max-min)*nb_segs );
	if ( seg<0 ) *below++;
	else if ( seg>nb_segs ) *above++;
	else distrib[ seg ]++;
    }
}	

int main(int argc, char **argv)
{
    if (argc!=5) return -1;
    float mu = atof(argv[1]),
	sigma = atof(argv[2]);
    unsigned nb_parts = atoi(argv[3]),
	nb_segs = atoi(argv[4]);
    float* parts,
	*rands_a,
	*rands_b;
    cudaMalloc(&parts, sizeof(float)*nb_parts);
    cudaMalloc(&rands_a, sizeof(float)*nb_parts);
    cudaMalloc(&rands_b, sizeof(float)*nb_parts);
    dim3 blockSize(512),
	gridSize(ceil((float) nb_parts/512));
    curandGenerator_t gen;
    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, time(NULL));
    curandGenerateUniform(gen, rands_a, nb_parts);
    curandGenerateUniform(gen, rands_b, nb_parts);
    trajs<<<gridSize, blockSize>>>(mu, sigma, parts, nb_parts, rands_a, rands_b);
    thrust::sort(thrust::device, parts, parts+nb_parts);
    unsigned* distrib,
	*above,
	*below;
    cudaMalloc(&distrib, sizeof(unsigned)*nb_segs);
    cudaMalloc(&below, sizeof(unsigned));
    cudaMalloc(&above, sizeof(unsigned));
    make_distrib<<<gridSize, blockSize>>>(parts,
					  nb_parts,
					  distrib,
					  nb_segs,
					  below,
					  above,
					  0, 1,
					  MIN(nb_segs/2, nb_parts/2));
    std::vector<unsigned> h_distrib (nb_segs);
    cudaMemcpy(h_distrib.data(), distrib, sizeof(unsigned)*nb_segs, cudaMemcpyDeviceToHost);
    // for (int i=0; i<nb_segs; i++)
    //     std::cout << (float) i/nb_segs << " " << h_distrib.at(i) << std::endl;
    return 0;
}