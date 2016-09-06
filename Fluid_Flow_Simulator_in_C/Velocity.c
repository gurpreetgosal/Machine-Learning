/*
   Update velocities u and v.
*/
#include <stdio.h>
#include <malloc.h>
#include <math.h>

int Velocity(int counter,int nx,int ny,double dt,double **F,double **G,double **p,
             double **u,double **v,double **phi,double del, double rho_L, double rho_A,double Re_L)
{
  int    i,j;
  double dx,dy;
  double rho_in, rho_out, rho_ave, theta_e, theta_n; 
  //printf("\n inside velocity\n");
  dx = 5.0/nx; dy = 2.0/ny;

  for (j=1; j<=ny; j++) {
    for (i=0; i<nx; i++) {
      //printf("\ninside u"); 
       if (phi[i][j] <=  0){
	 //printf("\n phi[i][j]<=0 and = %15.10f",phi[i][j]);
           rho_in  = rho_L;
           rho_out = rho_A;
          }
       else if (phi[i][j] >  0) {
	 //printf("\n phi[i][j] > 0and = %15.10f",phi[i][j]);
           rho_in  = rho_A;
           rho_out = rho_L;
	 }
      
       if (phi[i+1][j]>0 != phi[i][j]>0){
           theta_e = fabs(phi[i][j]/(phi[i][j] - phi[i+1][j]));
           rho_ave = rho_in*theta_e + rho_out*(1-theta_e);
	   // printf("\n rho_ave = %15.10f\ttheta_e = %15.10f",rho_ave,theta_e);
           //printf("\nfabs %15.10f",phi[i][j]/(phi[i][j] - phi[i+1][j]));
           u[i][j]=F[i][j]/del-(p[i+1][j]-p[i][j])*dt/(dx*rho_ave);              
           }
       else {
	 // printf("\n rho_in = %15.10f",rho_in);
	 u[i][j]=F[i][j]/del-(p[i+1][j]-p[i][j])*dt/(dx*rho_in);           
           }  
   
    }
  }
  for (i=1; i<=nx; i++) {
    for (j=1; j<ny; j++) {

      if (phi[i][j] <=  0){
	 //printf("\n phi[i][j]<=0 ");
           rho_in  = rho_L;
           rho_out = rho_A;
          }
       else if (phi[i][j] >  0) {
	 // printf("\n phi[i][j] > 0 ");
           rho_in  = rho_A;
           rho_out = rho_L;
	 }
      
       if (phi[i][j+1]>0 != phi[i][j]>0){
           theta_n = fabs(phi[i][j]/(phi[i][j] - phi[i][j+1]));
           rho_ave = rho_in*theta_n + rho_out*(1-theta_n);
           v[i][j]= (G[i][j]-(p[i][j+1]-p[i][j])*dt/(dy*rho_ave))/(del*del);  
           }
       else {
	   v[i][j]= (G[i][j]-(p[i][j+1]-p[i][j])*dt/(dy*rho_in))/(del*del);        
           }
      
    }
  }

  return 0;
}
