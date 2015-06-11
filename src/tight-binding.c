/**
 *
 * Tight-binding model for lead and Pt chain.
 * 
 * Gonzalo Aguirre <graguire@gmail.com>
 * Original idea: G. Chiappe
 */


/**
 * The input is a XYZ format file with the Pt atoms coordinates
 * The first stage an Hamiltonian will be created, atoms will be 
 *   related between each other if the square of distance is below
 *   certain threshold. The hopping between orbitals will be set 
 *   following the Slater-Koster's criteria. You can get the H 
 *   matrix so far.
 * The second stage calculate the eigen values vector.
 * The third stage is calculate the density of state, setting 
 *   min and max got from eigenvalues.
 */
#include "hamiltonian.h"
#include <omp.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

/*
#ifdef DEBUG
 #define DPRINT(ftm, arg...) fprintf(stderr, fmt, args...)
#else
 #define DPRINT(ftm, arg...) do{ } while(0) 
#endif
*/

#define sqr(x) ((x)*(x))
#define PI 3.1415926538
#define ORB 6
#define SPIN 2
#define ATMS 40


int	nl=36,
	na=4,
	nm=40*ORB;


/* numero cuantico de spin */
double ms_number[]={-1/2, 1/2}; 
/* numero cuantico ml para el orbital s: 0, y d: +2,...,-2, no cambiar el orden! */
double ml_number[]={0 ,2 ,1, 0, -1, -2}; 

/* cadena con el nombre de los orbitales con el orden original */
char *orbital_str[]={"s","xy","yz","z^2","xz","x^2-y^2"};


void usage(){
	printf("Usage:\n");
	printf("\t-d\tprint density of states\n");
	printf("\t-e\tprint eigenvalues\n");
	printf("\t-h\tprint hamiltonian matrix\n");
}


/* print square matrix */
void printMat(gsl_matrix * M, int N){
int i,j;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++)
			printf("%.3lf ", gsl_matrix_get(M,i,j));
		printf("\n");
	}


}

/* print square complex matrix */
void printComMat(gsl_matrix_complex * M, int N){
int i,j;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++)
			printf("%.3lf ", gsl_complex_abs(gsl_matrix_complex_get(M,i,j)));
		printf("\n");
	}


}

int main(int argc, char *argv[]){
int 	i,j,k,
	c,
	N;
int	dflag=0,
	eflag=0,
	hflag=0;

float 	w; /* frec */

//char *lvalue=NULL;

double 	**M, // XYZ coordinates
	dos,
	lambda=0;

	while((c = getopt (argc, argv, "dehl:")) != -1){
		switch (c){
			case 'd':
				dflag = 1;
				break;
			case 'e':
				eflag = 1;
				break;
			case 'h':
				hflag = 1;
				break;
			case 'l':
				lambda = atof(optarg);
				break;
		}
	}


	scanf("%d",&N);

	M = (double **) malloc (N*sizeof(double *)); // coordinate matrix

	// read coordinates (XYZ format file)
	for (int i=0; i<N; i++){
		char null[5]; // discard element
		double *tmp = (double *) malloc (3 * sizeof(double)); // 3 coordinates
		scanf("%s%lf%lf%lf", null, &tmp[0], &tmp[1], &tmp[2]);
		M[i] = tmp;
//		printf("- %.2f %.2f\n",M[i][0], M[i][1]); // DEBUG
	}

	/* M: coordinate matrix, N: number of atoms, l: spin-orbit parameter (set to 0 to tight-binding)*/
	gsl_matrix_complex * Hso = hamiltonian(M, N, lambda);
	/* print hamiltonial */
	if (hflag){
		printComMat(Hso,N*SPIN*ORB);
		return 0;
	}



	/* eigenvalues */
	gsl_matrix_complex * evec = gsl_matrix_complex_alloc(N*SPIN*ORB, N*SPIN*ORB);
	gsl_vector * eval = gsl_vector_alloc(N*SPIN*ORB);
	gsl_eigen_hermv_workspace * ws = gsl_eigen_hermv_alloc(N*SPIN*ORB);
	gsl_matrix_complex * A = Hso; // gsl_eigen_hermv() destroys Hso matrix, use a copy instead
	gsl_eigen_hermv (A, eval, evec, ws);
	gsl_eigen_hermv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_hermv_free(ws);

	if (eflag){
		for (int i=0; i<N*SPIN*ORB; i++)
			printf("%d %.4g \n", i, gsl_vector_get(eval,i));
		return 0;
	}


	/* calculate DoS */
	gsl_matrix_complex * G = gsl_matrix_complex_alloc(N*SPIN*ORB, N*SPIN*ORB); // Green
	gsl_matrix_complex_set_all(G, GSL_COMPLEX_ZERO); // inicializo

	double 	eval_min = gsl_vector_min (eval),
		eval_max = gsl_vector_max (eval);	

	for (w = eval_min; w < eval_max; w += 0.001){
		dos = 0;	
		#pragma omp parallel num_threads(4)
		{
		int tid = omp_get_thread_num();
		#pragma omp for private(i,k) reduction (+:dos)
		for (i=0; i<N*SPIN*ORB; i++)	
			for (k=0; k<N*SPIN*ORB; k++){
				gsl_complex h = gsl_matrix_complex_get (Hso, i, k);
				double l = gsl_vector_get (eval ,k);
				gsl_complex z = gsl_complex_rect(0,5e-3); /* parte imaginaria */
				gsl_complex num = gsl_complex_mul(h,gsl_complex_conjugate(h)); /* numerador */
				gsl_complex den = gsl_complex_add_real(z, w-l); /* denominador */
				gsl_complex g = gsl_complex_div(num,den);
				dos += GSL_IMAG(g);
			}
		if (dflag && tid==0)
			printf("%.3g %g \n", w, -dos/PI);
		}
	}
	
	gsl_matrix_complex_free(G);

	gsl_vector_free(eval);
	gsl_matrix_complex_free(evec);

	return 0;
}
