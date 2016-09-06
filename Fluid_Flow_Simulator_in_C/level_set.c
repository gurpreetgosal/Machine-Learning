//
//  level_set.c
//  gravity_currents
//
//  Created by Gurpreet Gosal on 2016-02-26.
//  Copyright © 2016 Gurpreet Gosal. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
define level_set */
int level_set(int nx,int ny,double dt,double **u,double **v,double **phi )
{    //printf("\ndone velocity\n");
    int    i,j,Nphi;
    double u_center, v_center, del_x, del_y,**new_phi;
    Nphi = (nx+2)*(ny+2);
    del_x = 5.0/nx; del_y = 2.0/ny;
    //Np = (nx+2)*(ny+2);
   //printf("\nentering level set\n %f\n%f\n",del_x,del_y);
     new_phi = (double **) malloc((nx+2)*sizeof(double *));
     new_phi[0] = (double *) malloc(Nphi*sizeof(double));
     for (i=1; i<nx+2; i++) new_phi[i]=new_phi[i-1]+(ny+2);
     for (i=0; i<Nphi; i++) new_phi[0][i]=0.0;
     /*
    for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
	  new_phi[i][j] = phi[i][j];
	}
	}*/
    
     
    for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            
            // Enforce Upwinding
            u_center = 0.5*(u[i-1][j] + u[i][j]);
            v_center = 0.5*(v[i][j-1] + v[i][j]);
            
            // FTBS in x and y
            if (u_center >= 0  && v_center >= 0 ) {
	      new_phi[i][j] = phi[i][j] - dt*( u_center*( phi[i][j] - phi[i-1][j] )/del_x + \
                                        v_center*( phi[i][j] - phi[i][j-1] )/del_y  );
            }
            // FTFS in x FTBS in y
            if (u_center < 0  && v_center >= 0 ) {
                new_phi[i][j] = phi[i][j] - dt*( u_center*( phi[i+1][j] - phi[i][j] )/del_x + \
                                            v_center*( phi[i][j] - phi[i][j-1] )/del_y  );
            }
            // FTFS in x FTFS in y
            if (u_center < 0  && v_center < 0 ) {
                new_phi[i][j] = phi[i][j] - dt*( u_center*( phi[i+1][j] - phi[i][j] )/del_x + \
                                            v_center*( phi[i][j+1] - phi[i][j] )/del_y  );
            }
            // FTBS in x FTFS in y
            if (u_center >= 0  && v_center < 0 ) {
                new_phi[i][j] = phi[i][j] - dt*( u_center*( phi[i][j] - phi[i-1][j] )/del_x + \
                                            v_center*( phi[i][j+1] - phi[i][j] )/del_y  );
            }
            
        }
    }
    copy_phi(nx,ny,phi,new_phi);
    phi_BC(nx,ny,phi);
    
    return 0;
}

int phi_BC(int nx, int ny, double **phi ){
    int i,j;
    /* set ghost values*/
    for (i=0; i<=nx+1; i++) {
        phi[i][0]= phi[i][1];
        phi[i][ny+1]=phi[i][ny];
    }
    for (j=0; j<=ny+1; j++) {
        phi[0][j]=phi[nx-1][j];
        phi[nx+1][j]=phi[2][j];
        phi[1][j]=phi[nx][j];
    }
    /*phi[0][0]=phi[nx+1][ny+1];
    phi[nx+1][0]=phi[0][ny+1];
    phi[0][ny+1]=phi[nx+1][0];
    phi[nx+1][ny+1]=phi[0][0];
    */
    return 0;
}

int copy_phi(int nx,int ny,double **phi,double **new_phi)
{
  int i,j;
for (j=0; j<=ny+1; j++) {
        for (i=0; i<=nx+1; i++) {
	  phi[i][j] = new_phi[i][j];
	}
 }
  return 0;
}
