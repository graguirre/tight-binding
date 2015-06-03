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
	/* operador LS */
//	gsl_complex z = gsl_complex_rect(0,1); /* i */
	gsl_matrix_complex * LS = gsl_matrix_complex_alloc(SPIN*ORB, SPIN*ORB); // operador LS
	gsl_matrix_complex_set_all(LS, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set(LS,1, 5,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(LS,2, 4,gsl_complex_rect(0,-0.5));
	gsl_matrix_complex_set(LS,4, 2,gsl_complex_rect(0, 0.5));
	gsl_matrix_complex_set(LS,5, 1,gsl_complex_rect(0, 1));

	gsl_matrix_complex_set(LS, 1, 8,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(LS, 1,10,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(LS, 2, 7,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(LS, 2, 9,gsl_complex_rect( 0         , sqrt(3)/2));
	gsl_matrix_complex_set(LS, 2,11,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(LS, 3, 8,gsl_complex_rect( 0         ,-sqrt(3)/2));
	gsl_matrix_complex_set(LS, 3,10,gsl_complex_rect(-sqrt(3)/2,  0));
	gsl_matrix_complex_set(LS, 4, 7,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(LS, 4, 9,gsl_complex_rect( sqrt(3)/2,  0));
	gsl_matrix_complex_set(LS, 4,11,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(LS, 5, 8,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(LS, 5,10,gsl_complex_rect( 0.5       , 0));

	gsl_matrix_complex_set(LS, 7, 2,gsl_complex_rect(-0.5       , 0));
	gsl_matrix_complex_set(LS, 7, 4,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(LS, 8, 1,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(LS, 8, 3,gsl_complex_rect( 0         , sqrt(3)/2));
	gsl_matrix_complex_set(LS, 8, 5,gsl_complex_rect( 0         , 0.5));
	gsl_matrix_complex_set(LS, 9, 2,gsl_complex_rect( 0         ,-sqrt(3)/2));
	gsl_matrix_complex_set(LS, 9, 4,gsl_complex_rect( sqrt(3)/2,  0));
	gsl_matrix_complex_set(LS,10, 1,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(LS,10, 3,gsl_complex_rect(-sqrt(3)/2,  0));
	gsl_matrix_complex_set(LS,10, 5,gsl_complex_rect( 0.5       , 0));
	gsl_matrix_complex_set(LS,11, 2,gsl_complex_rect( 0         ,-0.5));
	gsl_matrix_complex_set(LS,11, 4,gsl_complex_rect(-0.5       , 0));

	gsl_matrix_complex_set(LS, 7,11,gsl_complex_rect(0,-1));
	gsl_matrix_complex_set(LS, 8,10,gsl_complex_rect(0,-0.5));
	gsl_matrix_complex_set(LS,10, 8,gsl_complex_rect(0, 0.5));
	gsl_matrix_complex_set(LS,11, 7,gsl_complex_rect(0, 1));


	gsl_matrix_complex * Hso = gsl_matrix_complex_alloc(N*ORB*SPIN,N*ORB*SPIN); // spin-orbit Hamiltonian

	/* orbit order:  s , xy , yz , z^2 , xz , x^2-y^2 */
	/* spin order: up, down */
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for(int k0=0; k0<SPIN*ORB; k0++)		
			for(int k1=0; k1<SPIN*ORB; k1++){	
				//printf("%d %d %d %d\n",i,j,k0,k1);
				if (i==j){ // spin-orbit interaction
					gsl_complex z = gsl_matrix_complex_get(LS,k0,k1);
					gsl_matrix_complex_set(Hso, i*SPIN*ORB+k0, j*SPIN*ORB+k1, z); 
				}else { // tight-binding
					gsl_matrix_complex_set(Hso, i*SPIN*ORB+k0, j*SPIN*ORB+k1, gsl_complex_rect(0,1)); 
				}
			}

					
	
	/* print hamiltonial */
	if (hflag){
		printMat(H,N*ORB);
		return 0;
	}

	if (hflag){
		printComMat(Hso,N*SPIN*ORB);
		return 0;
	}


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
