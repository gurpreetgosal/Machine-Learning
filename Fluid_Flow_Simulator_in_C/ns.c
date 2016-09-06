#include <stdio.h>
#include <malloc.h>

/*
   Navier-Stokes solver:

    u_t = NS

*/

int main(int argc, char **args)
{
  int    nx,ny,i=1;
  double time,**u,**v,**p,**phi;
  printf("in ns.c about to enter Setup.c\n");
  /* Initialization */
  Setup("parm",&nx,&ny,&time,&u,&v,&p);
  printf("Setup finished");
  printf("nx=%d, ny=%d, time=%f\n",nx,ny,time);
  /* Solve the system of PDE */
  printf("entering NSSolve.c\n");
    
  NSSolve(nx,ny,time,u,v,p);

    
  /* Free everything */
  free(u[0]); free(u); free(v[0]); free(v);
  free(p[0]); free(p);

  return 0;
}

