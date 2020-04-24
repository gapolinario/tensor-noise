/*******************************************************
       File:    tensor-noise.c
       Name:    Gabriel B. Apolinario
      Usage:    ./tensor-noise.c
                describe parameters
    Compile:    -O3 -lm -lgsl -lgslcblas
     System:    ---
       Date:    24/04/2020
       Bugs:    No known bugs
Description:    
                Generate sequence of random tensors F_ij
                correlated as <F_ij(t) F_kl(t')> = G_ijkl \delta(t-t')

       Todo:    ---

*******************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// functions

int main(int argc, char **argv){

	int steps;
	double s;
	int i,ii; // loop integers

	steps = 10000;
	
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937); // Mersenne Twister
	gsl_rng_set(rng,100); // seed for random numbers

	gsl_vector * F = gsl_vector_alloc (9);
	gsl_matrix * B = gsl_matrix_alloc (9, 9); // B = Sqrt[G]
	gsl_vector * z = gsl_vector_alloc (9);
	gsl_vector * w = gsl_vector_alloc (9);

	char name[100];
    FILE *fout;
    sprintf(name, "data/random-tensors.dat");
    fout = fopen(name,"w");

	// set elements of B
	// they were calculated on Mathematica, RFD.nb
	// using LDL decomposition
	gsl_matrix_set(B,0,0,1.000000000000000);
	gsl_matrix_set(B,0,1,0);
	gsl_matrix_set(B,0,2,0);
	gsl_matrix_set(B,0,3,0);
	gsl_matrix_set(B,0,4,0);
	gsl_matrix_set(B,0,5,0);
	gsl_matrix_set(B,0,6,0);
	gsl_matrix_set(B,0,7,0);
	gsl_matrix_set(B,0,8,0);
	gsl_matrix_set(B,1,0,0);
	gsl_matrix_set(B,1,1,1.414213562373095);
	gsl_matrix_set(B,1,2,0);
	gsl_matrix_set(B,1,3,0);
	gsl_matrix_set(B,1,4,0);
	gsl_matrix_set(B,1,5,0);
	gsl_matrix_set(B,1,6,0);
	gsl_matrix_set(B,1,7,0);
	gsl_matrix_set(B,1,8,0);
	gsl_matrix_set(B,2,0,0);
	gsl_matrix_set(B,2,1,0);
	gsl_matrix_set(B,2,2,1.414213562373095);
	gsl_matrix_set(B,2,3,0);
	gsl_matrix_set(B,2,4,0);
	gsl_matrix_set(B,2,5,0);
	gsl_matrix_set(B,2,6,0);
	gsl_matrix_set(B,2,7,0);
	gsl_matrix_set(B,2,8,0);
	gsl_matrix_set(B,3,0,0);
	gsl_matrix_set(B,3,1,-0.3535533905932738);
	gsl_matrix_set(B,3,2,0);
	gsl_matrix_set(B,3,3,1.369306393762915);
	gsl_matrix_set(B,3,4,0);
	gsl_matrix_set(B,3,5,0);
	gsl_matrix_set(B,3,6,0);
	gsl_matrix_set(B,3,7,0);
	gsl_matrix_set(B,3,8,0);
	gsl_matrix_set(B,4,0,-0.5000000000000000);
	gsl_matrix_set(B,4,1,0);
	gsl_matrix_set(B,4,2,0);
	gsl_matrix_set(B,4,3,0);
	gsl_matrix_set(B,4,4,0.8660254037844386);
	gsl_matrix_set(B,4,5,0);
	gsl_matrix_set(B,4,6,0);
	gsl_matrix_set(B,4,7,0);
	gsl_matrix_set(B,4,8,0);
	gsl_matrix_set(B,5,0,0);
	gsl_matrix_set(B,5,1,0);
	gsl_matrix_set(B,5,2,0);
	gsl_matrix_set(B,5,3,0);
	gsl_matrix_set(B,5,4,0);
	gsl_matrix_set(B,5,5,1.414213562373095);
	gsl_matrix_set(B,5,6,0);
	gsl_matrix_set(B,5,7,0);
	gsl_matrix_set(B,5,8,0);
	gsl_matrix_set(B,6,0,0);
	gsl_matrix_set(B,6,1,0);
	gsl_matrix_set(B,6,2,-0.3535533905932738);
	gsl_matrix_set(B,6,3,0);
	gsl_matrix_set(B,6,4,0);
	gsl_matrix_set(B,6,5,0);
	gsl_matrix_set(B,6,6,1.369306393762915);
	gsl_matrix_set(B,6,7,0);
	gsl_matrix_set(B,6,8,0);
	gsl_matrix_set(B,7,0,0);
	gsl_matrix_set(B,7,1,0);
	gsl_matrix_set(B,7,2,0);
	gsl_matrix_set(B,7,3,0);
	gsl_matrix_set(B,7,4,0);
	gsl_matrix_set(B,7,5,-0.3535533905932738);
	gsl_matrix_set(B,7,6,0);
	gsl_matrix_set(B,7,7,1.369306393762915);
	gsl_matrix_set(B,7,8,0);
	gsl_matrix_set(B,8,0,-0.5000000000000000);
	gsl_matrix_set(B,8,1,0);
	gsl_matrix_set(B,8,2,0);
	gsl_matrix_set(B,8,3,0);
	gsl_matrix_set(B,8,4,-0.8660254037844386);
	gsl_matrix_set(B,8,5,0);
	gsl_matrix_set(B,8,6,0);
	gsl_matrix_set(B,8,7,0);
	gsl_matrix_set(B,8,8,0);

	for (ii = 0; ii < steps; ++ii)
	{

		// set random iid vector
		for (i = 0; i < 9; ++i)
	    	gsl_vector_set (z, i, gsl_ran_gaussian(rng, 1.));

	    // tensor product F_ij = B_ijkl z_kl
	    for (i = 0; i < 9; ++i)
	    {
			gsl_matrix_get_row(w, B, i);
			gsl_blas_ddot(w,z,&s);
			gsl_vector_set(F,i,s);
	    }

	    //gsl_vector_fprintf(stdout, F, "%.5g");
	    gsl_vector_fwrite(fout, F);

	}

    fclose(fout);

    gsl_vector_free(F);
    gsl_matrix_free(B);
    gsl_vector_free(z);
    gsl_vector_free(w);

	return 0;
}