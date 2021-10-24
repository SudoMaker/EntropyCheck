/*
    This file is part of EntropyCheck.

    Copyright (C) 2021 ReimuNotMoe <reimu@sudomaker.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the Apache License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    This file incorporates public domain works.
*/



#include "EntropyCheck.h"

#define _ISOC99_SOURCE
#include <math.h>

#define FALSE 0
#define TRUE  1

/*  RT_LOG2  --  Calculate log to the base 2  */
#define rt_log2 log2


/*

    Compute probability of measured Chi Square value.

    This code was developed by Gary Perlman of the Wang
    Institute (full citation below) and has been minimally
    modified for use in this program.

*/

#include <math.h>

/*HEADER
	Module:       z.c
	Purpose:      compute approximations to normal z distribution probabilities
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Copyright:    none
	Tabstops:     4
*/

#define	Z_MAX          6.0            /* maximum meaningful z value */

/*FUNCTION poz: probability of normal z value */
/*ALGORITHM
	Adapted from a polynomial approximation in:
		Ibbetson D, Algorithm 209
		Collected Algorithms of the CACM 1963 p. 616
	Note:
		This routine has six digit accuracy, so it is only useful for absolute
		z values < 6.  For z values >= to 6.0, poz() returns 0.0.
*/
static double        /*VAR returns cumulative probability from -oo to z */
poz(const double z)  /*VAR normal z value */
{
	double y, x, w;

	if (z == 0.0) {
		x = 0.0;
	} else {
		y = 0.5 * fabs(z);
		if (y >= (Z_MAX * 0.5)) {
			x = 1.0;
		} else if (y < 1.0) {
			w = y * y;
			x = ((((((((0.000124818987 * w
				    -0.001075204047) * w +0.005198775019) * w
				  -0.019198292004) * w +0.059054035642) * w
				-0.151968751364) * w +0.319152932694) * w
			      -0.531923007300) * w +0.797884560593) * y * 2.0;
		} else {
			y -= 2.0;
			x = (((((((((((((-0.000045255659 * y
					 +0.000152529290) * y -0.000019538132) * y
				       -0.000676904986) * y +0.001390604284) * y
				     -0.000794620820) * y -0.002034254874) * y
				   +0.006549791214) * y -0.010557625006) * y
				 +0.011630447319) * y -0.009279453341) * y
			       +0.005353579108) * y -0.002141268741) * y
			     +0.000535310849) * y +0.999936657524;
		}
	}
	return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

/*
	Module:       chisq.c
	Purpose:      compute approximations to chisquare distribution probabilities
	Contents:     pochisq()
	Uses:         poz() in z.c (Algorithm 209)
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Copyright:    none
	Tabstops:     4
*/

#define	LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define	BIGX           20.0         /* max value to represent exp (x) */
#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp(x))

/*FUNCTION pochisq: probability of chi sqaure value */
/*ALGORITHM Compute probability of chi square value.
	Adapted from:
		Hill, I. D. and Pike, M. C.  Algorithm 299
		Collected Algorithms for the CACM 1967 p. 243
	Updated for rounding errors based on remark in
		ACM TOMS June 1985, page 185
*/

static double pochisq(
	const double ax,    /* obtained chi-square value */
	const int df	    /* degrees of freedom */
)
{
	double x = ax;
	double a, y, s;
	double e, c, z;
	int even;	    	    /* true if df is an even number */

	if (x <= 0.0 || df < 1) {
		return 1.0;
	}

	a = 0.5 * x;
	even = (2 * (df / 2)) == df;
	if (df > 1) {
		y = ex(-a);
	}
	s = (even ? y : (2.0 * poz(-sqrt(x))));
	if (df > 2) {
		x = 0.5 * (df - 1.0);
		z = (even ? 1.0 : 0.5);
		if (a > BIGX) {
			e = (even ? 0.0 : LOG_SQRT_PI);
			c = log(a);
			while (z <= x) {
				e = log(z) + e;
				s += ex(c * z - a - e);
				z += 1.0;
			}
			return (s);
		} else {
			e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
			c = 0.0;
			while (z <= x) {
				e = e * (a / z);
				c = c + e;
				z += 1.0;
			}
			return (c * y + s);
		}
	} else {
		return s;
	}
}

/*  RT_INIT  --  Initialise random test counters.  */
void EntropyCheck_Init(EntropyCheckContext *ctx, bool binary_mode) {
	ctx->binary = binary_mode;	       /* Set binary / byte mode */

	/* Initialise for calculations */

	ctx->ent = 0.0;		       /* Clear entropy accumulator */
	ctx->chisq = 0.0;	       /* Clear Chi-Square */
	ctx->datasum = 0.0;	       /* Clear sum of bytes for arithmetic mean */

	ctx->mp = 0;		       /* Reset Monte Carlo accumulator pointer */
	ctx->mcount = 0; 	       /* Clear Monte Carlo tries */
	ctx->inmont = 0; 	       /* Clear Monte Carlo inside count */
	ctx->incirc = 65535.0 * 65535.0;/* In-circle distance for Monte Carlo */

	ctx->sccfirst = TRUE;	       /* Mark first time for serial correlation */
	ctx->scct1 = ctx->scct2 = ctx->scct3 = 0.0; /* Clear serial correlation terms */

	ctx->incirc = pow(pow(256.0, (double) (__ENTCHK_MONTEN / 2)) - 1, 2.0);

	for (int i = 0; i < 256; i++) {
		ctx->ccount[i] = 0;
	}

	ctx->totalc = 0;
}

/*  RT_ADD  --	Add one or more bytes to accumulation.	*/
void EntropyCheck_Feed(EntropyCheckContext *ctx, void *buf, size_t buf_len) {
	unsigned char *bp = buf;
	int oc, c, bean;

	while (bean = 0, (buf_len-- > 0)) {
		oc = *bp++;

		do {
			if (ctx->binary) {
				c = !!(oc & 0x80);
			} else {
				c = oc;
			}
			ctx->ccount[c]++;		  /* Update counter for this bin */
			ctx->totalc++;

			/* Update inside / outside circle counts for Monte Carlo
			   computation of PI */

			if (bean == 0) {
				ctx->monte[ctx->mp++] = oc;       /* Save character for Monte Carlo */
				if (ctx->mp >= __ENTCHK_MONTEN) {     /* Calculate every MONTEN character */
					int mj;

					ctx->mp = 0;
					ctx->mcount++;
					ctx->montex = ctx->montey = 0;
					for (mj = 0; mj < __ENTCHK_MONTEN / 2; mj++) {
						ctx->montex = (ctx->montex * 256.0) + ctx->monte[mj];
						ctx->montey = (ctx->montey * 256.0) + ctx->monte[(__ENTCHK_MONTEN / 2) + mj];
					}
					if ((ctx->montex * ctx->montex + ctx->montey *  ctx->montey) <= ctx->incirc) {
						ctx->inmont++;
					}
				}
			}

			/* Update calculation of serial correlation coefficient */

			ctx->sccun = c;
			if (ctx->sccfirst) {
				ctx->sccfirst = FALSE;
				ctx->scclast = 0;
				ctx->sccu0 = ctx->sccun;
			} else {
				ctx->scct1 = ctx->scct1 + ctx->scclast * ctx->sccun;
			}
			ctx->scct2 = ctx->scct2 + ctx->sccun;
			ctx->scct3 = ctx->scct3 + (ctx->sccun * ctx->sccun);
			ctx->scclast = ctx->sccun;
			oc <<= 1;
		} while (ctx->binary && (++bean < 8));
	}
}

/*  RT_END  --	Complete calculation and return results.  */
void EntropyCheck_Finalize(EntropyCheckContext *ctx, EntropyCheckResult *res) {
	int i;

	/* Complete calculation of serial correlation coefficient */

	ctx->scct1 = ctx->scct1 + ctx->scclast * ctx->sccu0;
	ctx->scct2 = ctx->scct2 * ctx->scct2;
	ctx->scc = ctx->totalc * ctx->scct3 - ctx->scct2;
	if (ctx->scc == 0.0) {
		ctx->scc = -100000;
	} else {
		ctx->scc = (ctx->totalc * ctx->scct1 - ctx->scct2) / ctx->scc;
	}

	/* Scan bins and calculate probability for each bin and
	   Chi-Square distribution.  The probability will be reused
	   in the entropy calculation below.  While we're at it,
	   we sum of all the data which will be used to compute the
	   mean. */

	ctx->cexp = ctx->totalc / (ctx->binary ? 2.0 : 256.0);  /* Expected count per bin */
	for (i = 0; i < (ctx->binary ? 2 : 256); i++) {
		double a = ctx->ccount[i] - ctx->cexp;;

		ctx->prob[i] = ((double) ctx->ccount[i]) / ctx->totalc;
		ctx->chisq += (a * a) / ctx->cexp;
		ctx->datasum += ((double) i) * ctx->ccount[i];
	}

	/* Calculate entropy */

	for (i = 0; i < (ctx->binary ? 2 : 256); i++) {
		if (ctx->prob[i] > 0.0) {
			ctx->ent += ctx->prob[i] * rt_log2(1 / ctx->prob[i]);
		}
	}

	/* Calculate Monte Carlo value for PI from percentage of hits
	   within the circle */
	ctx->montepi = 4.0 * (((double) ctx->inmont) / ctx->mcount);

	/* Return results through struct */
	res->totalc = ctx->totalc;
	res->ent = ctx->ent;
	res->chisq = ctx->chisq;
	res->mean = ctx->datasum / ctx->totalc;
	res->montepi = ctx->montepi;
	res->scc = ctx->scc;

	res->pochisq = pochisq(res->chisq, (ctx->binary ? 1 : 255));
}