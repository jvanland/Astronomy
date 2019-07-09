#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define PI 3.141592653589793  /* pi to unnecessary precision */

void fg(double mu,double *x,double *v,double dt) {
	double f,g,fd,gd;			/* Gauss's f, g, fdot and gdot */
	double r,vsq;
	double u;					/* r v cos(phi) */
	double a;					/* semi-major axis */
	double e;					/* eccentricity */
	double ec,es;				/* e cos(E), e sin(E) */
	double en;					/* mean motion */
	double nf;					/* orbital frequency */
	double dec;					/* delta E */
	double dm;					/* delta mean anomoly */
	double lo = -2*PI;
	double up = 2*PI;
	double w;					/* function to zero */
	double wp;					/* first derivative */
	double wpp;					/* second derivative */
	double wppp;				/* third derivative */
	double dx,s,c,sh;
	double next;
	double converge;			/* converge criterion */
	int iter,i,j;
	const double DOUBLE_EPS = 1.2e-16;
	const int MAX_ITER = 256;

	/* 
	 * Evaluate some orbital quantites.
	 */
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	u = x[0]*v[0] + x[1]*v[1] + x[2]*v[2];
	a = 1/(2/r-vsq/mu);
	en = sqrt(mu/(a*a*a));
	ec = 1-r/a;
	es = u/(en*a*a);
	e = sqrt(ec*ec + es*es);
	nf = en/(2*PI);
	j = (int)(nf*dt);
	dt -= j/nf;					/* reduce to single orbit */
	dm = en*dt - es;
	if ((es*cos(dm)+ec*sin(dm)) > 0) dec = dm + 0.85*e;
	else dec = dm - 0.85*e;
	dm = en*dt;					/* reset mean anomoly */
	converge = fabs(dm*DOUBLE_EPS);
	/*
	 * First double precision iterations for dec solution.
	 * This solves Kepler's equation in difference form:
	 * dM = dE - e cos(E_0) sin(dE) + e sin(E_0)(1 - cos(dE))
	 * Note that (1 - cos(dE)) == 2*sin(dE/2)*sin(dE/2).  The latter is
	 * probably more accurate for small dE.
	 */
	for(iter=1;iter<=MAX_ITER;iter++) {
		s = sin(dec);
		c = cos(dec);
		sh = sin(0.5*dec);
		w = dec - ec*s + es*2*sh*sh - dm;
		if (w > 0) up = dec;
		else lo = dec;
		wp = 1 - ec*c + es*s;
		wpp = ec*s + es*c;
		wppp = ec*c - es*s;
		dx = -w/wp;
		dx = -w/(wp + 0.5*dx*wpp);
		dx = -w/(wp + 0.5*dx*wpp + dx*dx*wppp/6);
		next = dec + dx;
		if (fabs(dx) <= converge) break;
		if (next > lo && next < up) dec = next;
		else dec = 0.5*(lo + up);
		if (dec==lo || dec==up) break;
		}
	if (iter>MAX_ITER) {
	  printf("dec soln, iter: %d dconverge: %g dx: %g dec: %g\n",
		 iter, converge, dx, dec);
	  printf("%g %g %g %g\n",dm, en, ec, es);
	}
	/*
	 * Update the orbit.
	 * See Danby, eq. 6.8.11 and 6.8.12, for expressions for f, g, fdot,
	 * and gdot.
	 * JS: changed f and gd expressions to take advantage of the half angle 
	 * formula given above.
	 */
	f = 1 - (a/r)*2*sh*sh;
	g = dt + (s-dec)/en;
	fd = -(a/(r*wp))*en*s;
	gd = 1 - 2*sh*sh/wp;
	for (i=0;i<3;i++) {
		s = f*x[i]+g*v[i];
		v[i] = fd*x[i]+gd*v[i];
		x[i] = s;
		}
	}
