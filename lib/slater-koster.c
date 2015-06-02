#include <math.h>
#include "slater-koster.h"

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
		{Esx2y2, Ezero, Ezero  , Ex2y2z2 , Ezero  , Ex2y2x2y2} 
};



