#include <math.h>

/*
   Determine step size dt.
*/

int TimeStep(int nx,int ny,double **u,double **v,double *dt,double Re_L,double Re_A,double del)
{
  int    i,j;
  double umax,vmax,hx,hy,a,b,c,min,aA,aL;
  double Re=1000.0,tau=0.5;


  umax=fabs(u[0][1]);
  for (j=1; j<=ny; j++) {
    for (i=0; i<nx+1; i++) {
      if (fabs(u[i][j])>umax) umax=fabs(u[i][j]);
    }
  }

  vmax=fabs(v[1][0]);
  for (i=1; i<=nx; i++) {
    for (j=0; j<ny+1; j++) {
      if (fabs(v[i][j])>vmax) vmax=fabs(v[i][j]);
    }
  }

  hx=5.0/nx; hy=2.0/ny;
  aL=0.5*Re_L*del/( 1.0/(hx*hx) + 1.0/(hy*hy) );
  aA=0.5*Re_A*del/(1.0/(hx*hx) + 1.0/(hy*hy));
  b=hx/umax;
  c=hy/vmax;

  if (aL<aA){
    min= aL;}
  else{
    min = aA;
  }

  // min=a;
  if (b<min) min=b;
  if (c<min) min=c;
  min *= tau;

  if (*dt>min) *dt=min;

  return 0;
}

