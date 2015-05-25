/**
 *
 * Tight-binding model for lead and Pt chain.
 * 
 * Gonzalo Aguirre <graguire@gmail.com>
 *
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

/* PARAM: sss, dds, ddp, ddd, sds */
/* ultimo indice indica tipo de enlace, s: sigma. p: pi, d: delta */
/* referencia: Handbook of the band structure of elemental solids. D.A. Papaconstantopoulos */
#define SSS -0.05907
#define DDS -0.03753
#define DDP  0.02719
#define DDD -0.00488
#define SDS -0.03709

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

/* devuelve la energia dado un coseno director (R_j - R_i = <l,m,n>) */
/* referencia PR 94, 1498. Slater, Koster */
double Ess(double l, double m, double n){ return SSS; }

double Esx2y2(double l, double m, double n){ return sqrt(3)/2*l*l*SDS; } /* E_s,x^2-y^2 */
double Esz2  (double l, double m, double n){ return (n*n - l*l/2)*SDS; }/* E_s,3z^2-r^2 */
	
double Exyxy    (double l, double m, double n){ return                         0*DDS +                   l*l*DDP +                   n*n*DDD; } /* E_xy,xy */
double Exyyz    (double l, double m, double n){ return                         0*DDS +                   l*n*DDP -                   l*n*DDD; } /* E_xy,yz */
double Exzx2y2  (double l, double m, double n){ return               3/2*n*l*l*l*DDS +         l*n*(1-2*l*l)*DDP -         l*n*(1-l*l/2)*DDD; } /* E_xz,x^2-y^2 */
double Exzz2    (double l, double m, double n){ return   sqrt(3)*n*l*(n*n-l*l/2)*DDS + sqrt(3)*l*n*(l*l-n*n)*DDP -     sqrt(3)/2*l*l*l*n*DDD; } /* E_xz,3z^2-r^2 */
double Ex2y2x2y2(double l, double m, double n){ return               3/4*l*l*l*l*DDS +           l*l*(1-l*l)*DDP +      n*n*+1/4*l*l*l*l*DDD; } /* E_x^2-y^2,x^2-y^2 */
double Ex2y2z2  (double l, double m, double n){ return sqrt(3)/2*l*l*(n*n-l*l/2)*DDS -       sqrt(3)*l*l*n*n*DDP + sqrt(3)/4*l*l*(1+n*n)*DDD; } /* E_x^2-y^2,3z^2-r^2 */
double Ez2z2    (double l, double m, double n){ return   (n*n-l*l/2)*(n*n-l*l/2)*DDS +             3*l*l*n*n*DDP +           3/4*l*l*l*l*DDD; } /* E_3z^2-r^2,3z^2-r^2 */

double Ezero     (double l, double m, double n){ return 0.0; }

/* matrix de funciones con la energia de Slater-Koster */
/* orden orbitales:  s , xy , yz , z^2 , xz , x^2-y^2 */
double (*func[][ORB])(double l, double m, double n) = { \
		{Ess   , Ezero, Ezero  , Esz2    , Ezero  , Esx2y2}, \
		{Ezero , Exyxy, Exyyz  , Ezero   , Ezero  , Ezero}, \
		{Ezero , Exyyz, Ezero  , Ezero   , Ezero  , Ezero}, \
		{Esz2  , Ezero, Ezero  , Ez2z2   , Exzz2  , Ex2y2z2}, \
		{Ezero , Ezero, Ezero  , Exzz2   , Ezero  , Ezero}, \
		{Esx2y2, Ezero, Ezero  , Ex2y2z2 , Ezero  , Ex2y2x2y2} };


void usage(){
	printf("Usage:\n");
	printf("\t-d\tprint density of states\n");
	printf("\t-e\tprint eigenvalues\n");
	printf("\t-h\tprint hamiltonian matrix\n");
}

void printMat(gsl_matrix * M){
int i,j;
	for (i=0; i<nm; i++){
		for (j=0; j<nm; j++)
			printf("%.5g ", gsl_matrix_get(M,i,j));
		printf("\n");
	}


}


gsl_matrix * hamiltonian(int n, double ** d){
gsl_matrix * m = NULL;

	return m;
}
double cos_dir(double n, double x, double y, double z){
	return n/sqrt(x*x+y*y+z*z);
}

int main(int argc, char *argv[]){
int 	c,
	i,
	j,
	k,
	N;
int	dflag=0,
	eflag=0,
	hflag=0;

float 	f[4],
	w; /* frec */

double 	M[44][2], // posiciones
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

	// read coordinates
	for (i=0; i<N; i++){
		for (j=0; j<4; j++) // read: _ X Y Z
			scanf("%f",&f[j]);
		M[i][0] = f[2]; // X
		// plane structure Y=0
		M[i][1] = f[1]; // Z
//		printf("- %.2f %.2f\n",M[i][0], M[i][1]); // DEBUG
	}


	gsl_matrix * H = gsl_matrix_alloc(nm,nm); // create  Hamiltonian matrix 
	gsl_matrix_set_zero(H); // inicialize
	
	// read leads
	for (i=0; i<nl; i++)	
		for (j=0; j<nl; j++){
			double l = M[i][0]-M[j][0]; 	// proyection X
			double m = 0;			// proyection Y
			double n = M[i][1]-M[j][1];	// proyection Z
			l = cos_dir(l, l, m, n);
			m = 0; //cos_dir(m, l, m, n);
			n = cos_dir(n, l, m, n);
			for (int k0=0; k0<ORB; k0++)
				for (int k1=0; k1<ORB; k1++) {
					if ( sqr(M[i][0]-M[j][0]) + sqr(M[i][1]-M[j][1]) < 3.24 && i!=j) { // square of the difference
						gsl_matrix_set (H, k0*(nl+na)+i, k1*(nl+na)+j, func[k0][k1] (l,m,n)); // set hopping
						// uncomment if debbuging
//						fprintf(stderr,"--------- %d %d -----------\n",i,j);
//						fprintf(stderr, "Orbitales (%.2f, %.2f): %s,%s h: %g \n", l, n, orbital_str[k0], orbital_str[k1], func[k0][k1] (l,m,n));
						
					}
				}
		}

	// pt chain
	for (i=0; i<na-1; i++){
		double l = M[(nl+na)-na+i][0]-M[(nl+na)-na+i+1][0];
		double m = 0;
		double n = M[(nl+na)-na+i][1]-M[(nl+na)-na+i+1][1];
		l = cos_dir(l, l, m, n);
		m = 0; //cos_dir(m, l, m, n);
		n = cos_dir(n, l, m, n);
		for (int k0=0; k0<ORB; k0++)
			for (int k1=0; k1<ORB; k1++){
				gsl_matrix_set (H, (k0+1)*(nl+na)-na+i, (k1+1)*(nl+na)-na+i+1, func[k0][k1] (l,m,n));
				gsl_matrix_set (H, (k0+1)*(nl+na)-na+i+1, (k1+1)*(nl+na)-na+i, func[k0][k1] (l,m,n)); // symmetry
//				fprintf(stderr, "Orbitales (%.2f, %.2f): %s-%s h: %g \n", l, n, orbital_str[k0], orbital_str[k1], func[k0][k1] (l,m,n));
			}
	}



	// lead -- chain hopping
	for (int k0=0; k0<ORB; k0++)
		for (int k1=0; k1<ORB; k1++){
			double l = M[17][0]-M[36][0];
			double m = 0;
			double n = M[17][1]-M[36][1];
			l = cos_dir(l, l, m, n);
			m = 0; //cos_dir(m, l, m, n);
			n = cos_dir(n, l, m, n);

			gsl_matrix_set (H, k0*(nl+na)+17, k1*(nl+na)+36, func[k0][k1] (l,m,n));
			gsl_matrix_set (H, k0*(nl+na)+36, k1*(nl+na)+17, func[k0][k1] (l,m,n));
//			fprintf(stderr,"--------- %d %d -----------\n",k0,k1);
//			fprintf(stderr, "Orbitales (%.2f, %.2f): %s-%s h: %g \n", l, n, orbital_str[k0], orbital_str[k1], func[k0][k1] (l,m,n));
			l = M[18][0]-M[39][0];
			m = 0;
			n = M[18][1]-M[39][1];
			l = cos_dir(l, l, m, n);
			m = 0; //cos_dir(m, l, m, n);
			n = cos_dir(n, l, m, n);

			gsl_matrix_set (H, k0*(nl+na)+18, k1*(nl+na)+39, func[k0][k1] (l,m,n));
			gsl_matrix_set (H, k0*(nl+na)+39, k1*(nl+na)+18, func[k0][k1] (l,m,n));
//			fprintf(stderr, "Orbitales (%.2f, %.2f): %s-%s h: %g \n", l, n, orbital_str[k0], orbital_str[k1], func[k0][k1] (l,m,n));
		}
	// energy (digonal values)
	for (i=0; i<(nl+na); i++)
		for (j=0; j<ORB; j++)
			gsl_matrix_set (H, j*(nl+na)+i, j*(nl+na)+i, PARAM_E[j]);
	

	/* ---------- interaccion SPIN-ORBIT (en base cartesiana) ----------- */
	/* Referencia: */
	/* operador Lz */
	gsl_complex z = gsl_complex_rect(0,1); /* i */
	gsl_matrix_complex * Lz = gsl_matrix_complex_alloc(6,6); // operador Lz
	gsl_matrix_complex * Lm = gsl_matrix_complex_alloc(6,6); // operador L-
	gsl_matrix_complex * Lp = gsl_matrix_complex_alloc(6,6); // operador L+
	gsl_matrix_complex_set_all(Lz, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set_all(Lm, GSL_COMPLEX_ZERO); // inicializo
	gsl_matrix_complex_set_all(Lp, GSL_COMPLEX_ZERO); // inicializo

/*	double Lz[6][6] = {{0, 0, 0, 0, 0, 0},{0, 2, 0, 0, 0, 0},{0, 0, -1, 0, 0, 0}, \
			   {0, 0, 0, 1, 0, 0},{0, 0, 0, 0, -2, 0},{0, 0, 0, 0, 0, 0}};
*/
	/* operador L- */
/*	double Lm[6][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 1, 0, 0},{0, 0, 0, 0, 1, 0}, \
			   {0, 0, 0, 0, 0, 1},{0, 0, 0, 0, 0, 0},{0, 0, 1, 0, 0, 0}};
*/
	/* operador L+ */
/*	double Lp[6][6] = {{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 1}, \
			   {0, 1, 0, 0, 0, 0},{0, 0, 1, 0, 0, 0},{0, 0, 0, 1, 0, 0}};
*/
	gsl_matrix * Hso = gsl_matrix_alloc(nm*SPIN,nm*SPIN); // spin-orbit Hamiltonian

	/* orden orbitales:  s , xy , yz , z^2 , xz , x^2-y^2 */
	for (i=0; i<ATMS; i++)
		for (j=0; j<ORB*SPIN; j++)
			gsl_matrix_set(Hso, i+j*ATMS, i+j*ATMS, 0.5*ml_number[j%6]); /* Lz */
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

					
	

	if (hflag){
		printMat(H);
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
		for (i=0; i<nm; i++)
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
		#pragma omp parallel num_threads(4)
		{
		int tid = omp_get_thread_num();
		#pragma omp for private(i,k) reduction (+:dos)
		for (i=0; i<nm; i++)	
			for (k=0; k<nm; k++){
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
