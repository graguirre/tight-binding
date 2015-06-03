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

#include "slater-koster.h"
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

/* Ess, Edd*/
double PARAM_E[]={.45649, .41348, .41348, .41348, .41348, .41348};

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
			printf("%.5g ", gsl_matrix_get(M,i,j));
		printf("\n");
	}


}


double cos_dir(double n, double x, double y, double z){
	return n/sqrt(x*x+y*y+z*z);
}

int main(int argc, char *argv[]){
int 	c,
	N;
int	dflag=0,
	eflag=0,
	hflag=0;

float 	w; /* frec */

double 	**M, // XYZ coordinates
	dos;

	while((c = getopt (argc, argv, "deh")) != -1){
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


	gsl_matrix * H = gsl_matrix_alloc(N*ORB, N*ORB); // create  Hamiltonian matrix (alloc: atoms x orbitals) 
	gsl_matrix_set_zero(H); // inicialize

	// Phys. Rev. 25,753
	const float pt_thresh=7.6518; // Angstrong. Square of distance to first neighbor (2.766A). 3.912A lattice constant.

	// read coordinates to create hamiltonian
	for (int i=0; i<N; i++)	
		for (int j=0; j<N; j++)
			if ( i!=j ){
				double l = M[i][0]-M[j][0]; 	// proyection X
				double m = M[i][1]-M[j][1];	// proyection Y
				double n = M[i][2]-M[j][2];	// proyection Z
				for (int k0=0; k0<ORB; k0++)
					for (int k1=0; k1<ORB; k1++) {
//						printf ("%lf\n", l*l + m*m + n*n ); // square of the difference
						if ( l*l + m*m + n*n < pt_thresh ) { // square of the difference
							l = cos_dir(l, l, m, n);
							m = cos_dir(m, l, m, n);
							n = cos_dir(n, l, m, n);
							gsl_matrix_set (H, i*ORB+k0, j*ORB+k1, func[k0][k1] (l,m,n)); // set hopping
							// uncomment if debbuging
//							fprintf(stderr,"--------- %d %d -----------\n",i,j);
//							fprintf(stderr, "Orbitales (%.2f, %.2f): %s,%s h: %g \n", l, n, orbital_str[k0], orbital_str[k1], func[k0][k1] (l,m,n));
						}
					}
			}

	// energy (digonal values)
	for (int i=0; i<N; i++)
		for (int j=0; j<ORB; j++)
			gsl_matrix_set (H, i*ORB+j, i*ORB+j, PARAM_E[j]);
	

	/* ---------- interaccion SPIN-ORBIT (en base cartesiana) ----------- */
	/* Referencia: */
	/* operador Lz */
//	gsl_complex z = gsl_complex_rect(0,1); /* i */
	gsl_matrix_complex * Lz = gsl_matrix_complex_alloc(6,6); // operador Lz
	gsl_matrix_complex_set_all(Lz, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set(Lz,1,5,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(Lz,2,4,gsl_complex_rect(0,-0.5));
	gsl_matrix_complex_set(Lz,4,2,gsl_complex_rect(0, 0.5));
	gsl_matrix_complex_set(Lz,5,1,gsl_complex_rect(0, 1));

	gsl_matrix_complex * Lm = gsl_matrix_complex_alloc(6,6); // operador L-
	gsl_matrix_complex_set_all(Lm, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set(Lm,1,2,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(Lm,1,4,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(Lm,2,1,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(Lm,2,3,gsl_complex_rect( 0         , sqrt(3)/2));
	gsl_matrix_complex_set(Lm,2,5,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(Lm,3,2,gsl_complex_rect( 0         ,-sqrt(3)/2));
	gsl_matrix_complex_set(Lm,3,4,gsl_complex_rect(-sqrt(3)/2,  0));
	gsl_matrix_complex_set(Lm,4,1,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(Lm,4,3,gsl_complex_rect( sqrt(3)/2,  0));
	gsl_matrix_complex_set(Lm,4,5,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(Lm,5,2,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(Lm,5,4,gsl_complex_rect( 0.5       , 0));

	gsl_matrix_complex * Lp = gsl_matrix_complex_alloc(6,6); // operador L+
	gsl_matrix_complex_set_all(Lp, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set(Lp,1,2,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(Lp,1,4,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(Lp,2,1,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(Lp,2,3,gsl_complex_rect( 0         , sqrt(3)/2));
	gsl_matrix_complex_set(Lp,2,5,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(Lp,3,2,gsl_complex_rect( 0         ,-sqrt(3)/2));
	gsl_matrix_complex_set(Lp,3,4,gsl_complex_rect( sqrt(3)/2,  0));
	gsl_matrix_complex_set(Lp,4,1,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(Lp,4,3,gsl_complex_rect(-sqrt(3)/2,  0));
	gsl_matrix_complex_set(Lp,4,5,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(Lp,5,2,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(Lp,5,4,gsl_complex_rect(-0.5       , 0));

	gsl_matrix_complex * Hso = gsl_matrix_complex_alloc(N*ORB*SPIN,N*ORB*SPIN); // spin-orbit Hamiltonian

	/* orden orbitales:  s , xy , yz , z^2 , xz , x^2-y^2 */
	for (int i=0; i<N; i++)
		for (int j=0; j<ORB*SPIN; j++)
			gsl_matrix_complex_set(Hso, i+j*ATMS, i+j*ATMS, gsl_complex_rect(0,1)); /* Lz */
/*			
	for (i=0; i<ORB; i++)
		for (j=0; j<ORB; j++)
			if (Lm[i][j])
				for (k=0; k<ATMS; k++)
					gsl_matrix_set(Hso, ATMS*ORB+i*ATMS+k, i*ATMS+k, sqrt(2)); // L- 
			if (Lp[i][j])
				for (k=0; k<ATMS; k++)
					gsl_matrix_set(Hso, i*ATMS+k, ATMS*ORB+i*ATMS+k, sqrt(2)); // L+ 
*/

					
	
	/* print hamiltonial */
	if (hflag){
		printMat(H,N*ORB);
		return 0;
	}

	/* calculate eigenvalues */
	gsl_matrix * evec = gsl_matrix_alloc(nm,nm);
	gsl_vector * eval = gsl_vector_alloc(nm);
	gsl_eigen_symmv_workspace * ws = gsl_eigen_symmv_alloc(nm);
	gsl_eigen_symmv (H, eval, evec, ws);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free(ws);

	if (eflag){
		for (int i=0; i<nm; i++)
			printf("%d %.4g \n", i, gsl_vector_get(eval,i));
		return 0;
	}



	/* calculate DoS */
	gsl_matrix_complex * G = gsl_matrix_complex_alloc(nm,nm); // Green
	gsl_matrix_complex_set_all(G, GSL_COMPLEX_ZERO); // inicializo

	double eval_min = gsl_vector_min (eval),
	eval_max = gsl_vector_max (eval);	

	for (w = eval_min; w < eval_max; w += 0.001){
		dos = 0;	
//		#pragma omp parallel num_threads(4)
		{
		int tid = omp_get_thread_num();
//		#pragma omp for private(i,k) reduction (+:dos)
		for (int i=0; i<nm; i++)	
			for (int k=0; k<nm; k++){
				double h = gsl_matrix_get (H, i, k);
				double l = gsl_vector_get (eval ,k);
				gsl_complex z = gsl_complex_rect(0,1e-4); /* parte imaginaria */
				gsl_complex num = gsl_complex_rect(h*h,0); /* numerador */
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
	gsl_matrix_free(evec);

	gsl_matrix_free(H);

	return 0;
}
