/**
 * Return Slater-Koster parameters
 */

#ifndef _SLATER_KOSTER_H_
#define _SLATER_KOSTER_H_

#define ORB 6

/* PARAM: sss, dds, ddp, ddd, sds */
/* ultimo indice indica tipo de enlace, s: sigma. p: pi, d: delta */
/* referencia: Handbook of the band structure of elemental solids. D.A. Papaconstantopoulos */
#define SSS -0.05907
#define DDS -0.03753
#define DDP  0.02719
#define DDD -0.00488
#define SDS -0.03709



/* bindinb between orbitals */
double Ezero    (double , double , double );
double Ess      (double , double , double );
double Esx2y2   (double , double , double ); 
double Esz2     (double , double , double ); 
double Exyxy    (double , double , double ); 
double Exyyz    (double , double , double );
double Exzx2y2  (double , double , double );
double Exzz2    (double , double , double );
double Ex2y2x2y2(double , double , double );
double Ex2y2z2  (double , double , double );
double Ez2z2    (double , double , double );


/* matrix de funciones con la energia de Slater-Koster */
/* orden orbitales:  s , xy , yz , z^2 , xz , x^2-y^2 */
double (*func[6][ORB])(double , double , double );


#endif // _SLATER_KOSTER_H_
