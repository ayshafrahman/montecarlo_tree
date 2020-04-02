#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/***************************************************************
* Program to calculate the branches of a tree.
*
*	Units:
*		Length = m
*		Time   = sec
*
***************************************************************/

int main(int argc, char *argv[])
{
	double tpi = 6.283185308;
	double pi = 3.141592654;
	int i,j,k=0,m,n,p,q,t=0,s,l=0,Nlocal, Lsteps, rsteps;
	int MAXPTS=64000,MAXLOCAL=1000;
	FILE *outf;
	double stepsize=0.1, L0=10.0, r1, r2,r3,r4,t0=1;
	double A00, A01, A02, A10, A11, A12, A20, A21, A22;
	printf("here \n"); /*print*/
	double x0[2047], y0[2047], z0[2047], L[2047], thickness[2047];
	double cost[2047],sint[2047],phi[2047],cosp[2047],sinp[2047];
	double xg[MAXPTS], yg[MAXPTS], zg[MAXPTS];
	double xl[MAXLOCAL], yl[MAXLOCAL], zl[MAXLOCAL];
	char outfile[128];

	/* initialize */
	srand(time((time_t*)NULL));
	for(i=0; i<2048; i++)
	{
		L[i] = 0.0;
		thickness[i] = 0.0;
		x0[i] = 0.0;
		y0[i] = 0.0;
		z0[i] = 0.0;
	}
	for(i=0; i<MAXPTS; i++)
	{
		xg[i] = -1.0;
		yg[i] = -1.0;
		zg[i] = -1.0;
	}
	L[0] = L0; /* this could be user-specified */
	thickness[0] = t0;
	k=1;
	while(k < 2048)
	{
	/* generate a new branch */

	/* branch length */
		r1 = 1.0*rand()/RAND_MAX;
		L[k]=fabs(L[k-1]*r1);
		L[k+1] = L[k];
		r2 = 1.0*rand()/RAND_MAX;
		thickness[k]=fabs(thickness[k-1]*r2);
		thickness[k+1] = thickness[k-1]*(1-r2);
		Lsteps = L[k]/stepsize+1;
		rsteps = thickness[k]/stepsize+1;
			for(m=0;m<rsteps;m++)
			{
			for(p=0;p<rsteps;p++)
			{
			for(q=0;q<Lsteps;q++)
			{
			xl[l] = -thickness[k]+m*stepsize;
			yl[l] = -thickness[k]+p*stepsize;
			zl[l] = q*stepsize;
			l++;
			}}}
		Nlocal = Lsteps*rsteps*rsteps +1;

		/* compute branch direction */
			r3 = 1.0*rand()/RAND_MAX;
			/* direction cosine of z coordinate 0-pi/2 */
			cost[k]=r3;
			sint[k]=sqrt(1.0-cost[k]*cost[k]);
			cost[k+1] = cost[k];
			sint[k+1] = sint[k];
			/* one branch is the mirror image of the other -- off by pi */
			r4 = 1.0*rand()/RAND_MAX;
			phi[k]=tpi*r4;
			phi[k+1] = phi[k]-pi;
			cosp[k] = cos(phi[k]);
			sinp[k] = sin(phi[k]);

		/*  Rotation matrix */
		A00 = cosp[k]*cost[k];
		A01 = -sinp[k];
		A02 = cosp[k]*sint[k];
		A10 = sinp[k]*cost[k];
		A11 = cosp[k];
		A12 = sinp[k]*sint[k];
		A20 = -sint[k];
		A21 = 0.0;
		A22 = cost[k];
		for(l=0;l<Nlocal;l++)
		{
		xg[t] = A00*xl[l] + A01*yl[l] + A02*zl[l] + x0[k];
		yg[t] = A10*xl[l] + A11*yl[l] + A12*zl[l] + y0[k];
		zg[t] = A20*xl[l] + A21*yl[l] + A22*zl[l] + z0[k];
		t++;
		}

		x0[k+1] = A02*L[k] + x0[k];
		y0[k+1] = A12*L[k] + y0[k];
		z0[k+1] = A22*L[k] + z0[k];
		k++;
		}
/* check if finished */
    printf("Finished.\n");
}
