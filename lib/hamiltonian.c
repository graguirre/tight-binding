#include "hamiltonian.h"

double cos_dir(double n, double x, double y, double z){
	return n/sqrt(x*x+y*y+z*z);
}

gsl_matrix_complex * hamiltonian(double **M, int N, double lambda){


	/* LS operator (s + d orbitals)*/
	/* Ref.:  Jaffe, M. D., Singh, J. (1987/05)."Inclusion of spin-orbit coupling into tight binding bandstructure calculations for bulk and superlattice semiconductors." Solid State Communications 62(6): 399-402.  */
	gsl_matrix_complex * LS = gsl_matrix_complex_alloc(SPIN*ORB, SPIN*ORB); 
	gsl_matrix_complex_set_all(LS, GSL_COMPLEX_ZERO); // initialize
	gsl_matrix_complex_set(LS,1, 5,gsl_complex_rect(lambda*0,lambda*-1));
	gsl_matrix_complex_set(LS,2, 4,gsl_complex_rect(lambda*0,lambda*-0.5));
	gsl_matrix_complex_set(LS,4, 2,gsl_complex_rect(lambda*0,lambda* 0.5));
	gsl_matrix_complex_set(LS,5, 1,gsl_complex_rect(lambda*0,lambda* 1));

	gsl_matrix_complex_set(LS, 1, 8,gsl_complex_rect(lambda* 0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS, 1,10,gsl_complex_rect(lambda* 0         ,lambda* 0.5));
	gsl_matrix_complex_set(LS, 2, 7,gsl_complex_rect(lambda*-0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS, 2, 9,gsl_complex_rect(lambda* 0         ,lambda* sqrt(3)/2));
	gsl_matrix_complex_set(LS, 2,11,gsl_complex_rect(lambda* 0         ,lambda* 0.5));
	gsl_matrix_complex_set(LS, 3, 8,gsl_complex_rect(lambda* 0         ,lambda*-sqrt(3)/2));
	gsl_matrix_complex_set(LS, 3,10,gsl_complex_rect(lambda*-sqrt(3)/2 ,lambda*  0));
	gsl_matrix_complex_set(LS, 4, 7,gsl_complex_rect(lambda* 0         ,lambda*-0.5));
	gsl_matrix_complex_set(LS, 4, 9,gsl_complex_rect(lambda* sqrt(3)/2 ,lambda*  0));
	gsl_matrix_complex_set(LS, 4,11,gsl_complex_rect(lambda*-0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS, 5, 8,gsl_complex_rect(lambda* 0         ,lambda*-0.5));
	gsl_matrix_complex_set(LS, 5,10,gsl_complex_rect(lambda* 0.5       ,lambda* 0));

	gsl_matrix_complex_set(LS, 7, 2,gsl_complex_rect(lambda*-0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS, 7, 4,gsl_complex_rect(lambda* 0         ,lambda* 0.5));
	gsl_matrix_complex_set(LS, 8, 1,gsl_complex_rect(lambda* 0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS, 8, 3,gsl_complex_rect(lambda* 0         ,lambda* sqrt(3)/2));
	gsl_matrix_complex_set(LS, 8, 5,gsl_complex_rect(lambda* 0         ,lambda* 0.5));
	gsl_matrix_complex_set(LS, 9, 2,gsl_complex_rect(lambda* 0         ,lambda*-sqrt(3)/2));
	gsl_matrix_complex_set(LS, 9, 4,gsl_complex_rect(lambda* sqrt(3)/2 ,lambda*  0));
	gsl_matrix_complex_set(LS,10, 1,gsl_complex_rect(lambda* 0         ,lambda*-0.5));
	gsl_matrix_complex_set(LS,10, 3,gsl_complex_rect(lambda*-sqrt(3)/2 ,lambda*  0));
	gsl_matrix_complex_set(LS,10, 5,gsl_complex_rect(lambda* 0.5       ,lambda* 0));
	gsl_matrix_complex_set(LS,11, 2,gsl_complex_rect(lambda* 0         ,lambda*-0.5));
	gsl_matrix_complex_set(LS,11, 4,gsl_complex_rect(lambda*-0.5       ,lambda* 0));

	gsl_matrix_complex_set(LS, 7,11,gsl_complex_rect(lambda*0,lambda* 1));
	gsl_matrix_complex_set(LS, 8,10,gsl_complex_rect(lambda*0,lambda* 0.5));
	gsl_matrix_complex_set(LS,10, 8,gsl_complex_rect(lambda*0,lambda*-0.5));
	gsl_matrix_complex_set(LS,11, 7,gsl_complex_rect(lambda*0,lambda*-1));


	gsl_matrix_complex * Hso = gsl_matrix_complex_alloc(N*ORB*SPIN,N*ORB*SPIN); // spin-orbit Hamiltonian	
	// Phys. Rev. 25,753
	const float pt_thresh=7.8; // Angstrong. Square of distance to first neighbor (2.766A). 3.912A lattice constant.
	/* orbit order:  s , xy , yz , z^2 , xz , x^2-y^2 */
	/* spin order: up, down */
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			for(int k0=0; k0<SPIN*ORB; k0++)		
			for(int k1=0; k1<SPIN*ORB; k1++){	
//				printf("%d %d %d %d\n",i,j,k0,k1);
				if (i==j) { // spin-orbit interaction
					gsl_complex z = gsl_matrix_complex_get(LS,k0,k1);
					gsl_matrix_complex_set(Hso, i*SPIN*ORB+k0, j*SPIN*ORB+k1, z); 
				} else if ((k0 < ORB && k1 < ORB) || (k0 >= ORB && k1 >= ORB)){ // tight-binding (affect same spin)
					double l = M[i][0]-M[j][0]; 	// proyection X
					double m = M[i][1]-M[j][1];	// proyection Y
					double n = M[i][2]-M[j][2];	// proyection Z

					if ( l*l + m*m + n*n < pt_thresh ) { // square of the difference
//						printf("%.2f %.2f %.2f\n", l,m,n);
						l = cos_dir(l, l, m, n);
						m = cos_dir(m, l, m, n);
						n = cos_dir(n, l, m, n);
						double d = func[k0%6][k1%6] (l,m,n); // slater-koster 6x6 orbital matrix
						gsl_matrix_complex_set (Hso, i*SPIN*ORB+k0, j*SPIN*ORB+k1, gsl_complex_rect(d,0)); // set hopping
					}
				}
			}

	/* Ess, Edd*/
	double PARAM_E[]={.45649, .41348, .41348, .41348, .41348, .41348};
        // energy (digonal values)
        for (int i=0; i<N*SPIN*ORB; i++)
                gsl_matrix_complex_set (Hso, i, i, gsl_complex_rect(PARAM_E[i%6],0));


	return Hso;

}
