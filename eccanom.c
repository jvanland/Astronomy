/*
 ** This function solves Kepler's equation, M = E - e*sin(E), for E.
 ** It is mainly called by conversion from the Delaunay elements to 
 ** cartesian coordinates.
 */

#include <math.h>
#include <assert.h>
#include <stdio.h>

#define MAX_KEPLER_ITTR 50
#define PI 3.14159265358979323
#define THRESH 1.e-14
#define CUBE_ROOT( X)  (exp( log( X) / 3.))

double arcsinh( const double z)
{
   return( log( z + sqrt( z * z + 1.)));
}

double dEccAnom(double mean_anom,double ecc)
{

#ifdef Old_EccAnom   /*  Original Code to Solve Kepler's Equation  */
 /*
 ** Original code from the Saha and Tremaine long term planetary integrator. 
 **
 ** Joachim Stadel, Jan. 11, 1995
 */
        double w,w1,w2,w3,twopi,dummy;
#ifdef __i386__
	long double dE,E,lb,ub;
#else
	double dE,E,lb,ub;
#endif
	int i;

	twopi = 2.0*M_PI;
	mean_anom -= ((int)(mean_anom/twopi))*twopi;
	if (sin(mean_anom) > 0.0) E = mean_anom - 0.85*ecc;  /* sin(M) is -ve */
	else E = mean_anom + 0.85*ecc;  /* sin(M) is +ve */
	lb = -twopi;
	ub = twopi;
	for (i=0;i<MAX_KEPLER_ITTR;++i) {
		w2 = ecc*sin(E);
		w3 = ecc*cos(E);
		w = mean_anom + w2 - E;
		w1 = 1.0 - w3;
		dE = w/w1;
		dE = w/(w1 + 0.5*dE*w2);
		dE = w/(w1 + dE*(0.5*w2 + dE*w3/6.0));
		dummy = dE;
		if (dE < 0.0) ub = E;
		else lb = E;
		dE += E;
		if (E == dE) return(E);
		if (dE > lb && dE < ub) E = dE;
		else E = 0.5*(lb + ub);
		if (E == lb || E == ub) return(E);
		}
	/*
	 ** Kepler solution has failed.
	 */
	printf("M = %g,e = %g, E = %.16g, dE = %.16g\n",mean_anom,ecc,dE,dummy); 
	assert(0);
	return(0.0);

#else  /* New code for solving Kepler's Equation   */

   double curr, err;
   int is_negative = 0, n_iter = 0;

   if( !mean_anom) {
      return( 0.);
   }

   if( ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
     {
      curr = atan2( sin( mean_anom), cos( mean_anom) - ecc);
      err = curr - ecc * sin( curr) - mean_anom;
            /* one correction step,  and we're done */
      while( fabs( err) > THRESH)
	{
	  n_iter++;
	  if (n_iter >= MAX_KEPLER_ITTR){
	    printf("M = %g,e = %g, E = %.16g, dE = %.16g\n",mean_anom,ecc,curr,err);
	    assert(0);
	  }
	  curr -= err / (1. - ecc * cos( curr));
	  err = curr - ecc * sin( curr) - mean_anom;
	}
      return curr;
      }

   if( mean_anom < 0.)
      {
      mean_anom = -mean_anom;
      is_negative = 1;
      }

   curr = mean_anom;
   if( (ecc > 0.8 && mean_anom < PI / 3.0) || ecc > 1.0)    /* up to 60 degrees */
      {
      double trial = mean_anom / fabs( 1. - ecc);

      if( trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
         {
         if( mean_anom < PI)
            trial = CUBE_ROOT( 6. * mean_anom);
         else        /* hyperbolic w/ 5th & higher-order terms predominant */
            trial = arcsinh( mean_anom / ecc);
         }
      curr = trial;
      }

   if( ecc < 1.)
      {
      err = curr - ecc * sin( curr) - mean_anom;
      while( fabs( err) > THRESH)
         {
         n_iter++;
	 if (n_iter >= MAX_KEPLER_ITTR){
	   printf("M = %g,e = %g, E = %.16g, dE = %.16g\n",mean_anom,ecc,curr,err);
	   assert(0);
	 }
         curr -= err / (1. - ecc * cos( curr));
         err = curr - ecc * sin( curr) - mean_anom;
         }
      }
   else
      {
      err = ecc * sinh( curr) - curr - mean_anom;
      while( fabs( err) > THRESH)
         {
         n_iter++;
	 if (n_iter >= MAX_KEPLER_ITTR){
	   printf("M = %g,e = %g, E = %.16g, dE = %.16g\n",mean_anom,ecc,curr,err);
	   assert(0);
	 }
         curr -= err / (ecc * cosh( curr) - 1.);
         err = ecc * sinh( curr) - curr - mean_anom;
         }
      }

   return( is_negative ? -curr : curr);

#endif
	}



