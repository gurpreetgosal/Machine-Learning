#include<stdio.h>
/*
   Setup boundary conditions for u and v.
*/

int BC(int nx,int ny,double **u,double **v)
{
  int    i,j;

  printf("\n inside BC.c %d",nx);
  // printf("\n 2 BC");
  /* u Top/Bottom: free slip boundary condition */
  for (i=0; i<nx; i++) {
    u[i][0]=-u[i][1];
    u[i][ny+1]=u[i][ny];
    }

  /* u Periodic lateral boundary condition */
  for (j=0; j<=ny+1; j++) {
    // u[0][j]=u[nx-1][j];
    u[nx][j]=u[0][j];
    }

  /* No/free slip boundary condition for v*/
  for (i=0; i<=nx+1; i++) {
    v[i][0]=0.0;
    v[i][ny]=v[i][ny-1];
    //v[i][ny]=0.0;
    }
 /* v Periodic lateral boundary condition */
 for (j=0; j<ny+1; j++) {
    v[0][j]=v[nx-1][j];
    v[1][j]=v[nx][j];
    v[nx+1][j]=v[2][j]; 
    }
 //printf("\n 3 BC");
  return 0;
}
