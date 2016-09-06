//
//  RHS_final.c
//  gravity_currents
//
//  Created by Gurpreet Gosal on 2016-06-16.
//  Copyright © 2016 Gurpreet Gosal. All rights reserved.
//

#include <math.h>
#include <stdio.h>

/*
 Setup right-hand sides for u, v, p.
 */

int RHS(int nx,int ny,double dt,double **u,double **v,double **F,
        double **G,double **H,double **phi,double del,double rho_L,
        double rho_A, double mu_L,double mu_A,double beta)
{
    int    i,j;
    double d2udx2,d2udy2,du2dx,duvdy,gx;
    double d2vdx2,d2vdy2,duvdx,dv2dy,gy;
    double dx,dy,u1,u2,u3,u4,v1,v2,v3,v4;
    double phi_current_ux,phi_prev_ux,phi_next_ux,phi_current_uy,phi_prev_uy,phi_next_uy;
    double phi_current_vx,phi_prev_vx,phi_next_vx,phi_current_vy,phi_prev_vy,phi_next_vy;
    double theta_e,theta_w,theta_n,theta_s,mu_hat;
    double L = 1.0;
    //double Re_L = rho_L*1.0*1.0/mu_L , Re_A = rho_A*1.0*1.0/mu_A;
    double Re_L = 1.0 , Re_A = 1.0/15.0;
    double gamma=0.9;
    
    dx = 5.0/nx; dy = 2.0/ny;
    
    /* Compute F */
    for (j=1; j<=ny; j++) {
      //F[0][j]=del*u[0][j];
      //F[nx][j]=del*u[nx][j];
        /* i=0; i=nx */
        for (i=0; i<=nx; i++){
        // Computing exceptional cases i.e. F at left and right boundaries
        
            if (i==0) {
                
                
                u1=(u[i][j]+u[i+1][j])/2;
                u2=(u[nx-1][j]+u[i][j])/2;
                du2dx=(u1*u1-u2*u2)/dx;
                u1=fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])/4;
                u2=fabs(u[nx-1][j]+u[i][j])*(u[nx-1][j]-u[i][j])/4;
                du2dx += gamma*(u1-u2)/dx;
                //printf("\ndu2dx[0][1]=%15.10f",du2dx);
                v1=(v[i][j]+v[i+1][j])/2;
                v2=(v[i][j-1]+v[i+1][j-1])/2;
                u1=(u[i][j]+u[i][j+1])/2;
                u2=(u[i][j-1]+u[i][j])/2;
                u3=(u[i][j]-u[i][j+1])/2;
                u4=(u[i][j-1]-u[i][j])/2;
                duvdy=(v1*u1-v2*u2)/dy+gamma*(fabs(v1)*u3-fabs(v2)*u4)/dy;
                //printf("\n u[0][1]=%15.10f\tu[0][1]=%15.10f\tu[0][1]=%15.10f\tu[0][1]=%15.10f\t")
                //printf("\nduvdy[0][1]=%15.10f",duvdy);
                //gx =  3/Re_L;
                
                //determine phi for the
                phi_current_ux = (phi[i][j] + phi[i+1][j])/2.0;
                phi_next_ux    = (phi[i+1][j] + phi[i+2][j])/2.0;
                //phi_prev_ux    = (phi[i-1][j] + phi[i][j])/2.0;
                
                phi_current_uy = (phi[i][j] + phi[i+1][j])/2.0;
                phi_next_uy    = (phi[i][j+1] + phi[i+1][j+1])/2.0;
                phi_prev_uy    = (phi[i][j-1] + phi[i+1][j-1])/2.0;
                
                
                // determine d2udx2
                if (phi_current_ux  <=  0){
                    // in liquid
                    if ( phi_next_ux>0){
                        theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                        mu_hat = mu_A*theta_e + mu_L*(1-theta_e);
                        d2udx2 = (mu_A/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[nx-1][j])/(dx*dx);
                    }
                    /* else if  (phi_prev_ux>0){
                     theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                     mu_hat = mu_A*theta_w + mu_L*(1-theta_w);
                     d2udx2 = (u[i+1][j]-u[i][j])/(dx*dx)-(mu_A/mu_hat)*(u[i][j]-u[i-1][j])/(dx*dx);
                     }*/
                    else {
                        d2udx2=(u[nx-1][j]-2*u[i][j]+u[i+1][j])/(dx*dx);
                    }
                    
                }
                else if  (phi_current_ux >  0){
                    // in air
                    if ( phi_next_ux <= 0){
                        theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                        mu_hat = mu_L*theta_e + mu_A*(1-theta_e);
                        d2udx2 = (mu_L/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[nx-1][j])/(dx*dx);
                    }
                    /* else if  (phi_prev_ux <= 0){
                     theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                     mu_hat = mu_L*theta_w + mu_A*(1-theta_w);
                     d2udx2 = (u[i+1][j]-u[i][j])/(dx*dx)-(mu_L/mu_hat)*(u[i][j]-u[nx][j])/(dx*dx);
                     }*/
                    else {
                        d2udx2=(u[nx-1][j]-2*u[i][j]+u[i+1][j])/(dx*dx);
                    }
                    
                }
               // printf("\n\nd2udx2[0][1]=%15.10f",d2udx2);
                // determine d2vdy2
                if (phi_current_uy  <=  0){
                    // in liquid
                    if ( phi_next_uy>0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_A*theta_n + mu_L*(1-theta_n);
                        d2udy2 = (mu_A/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy>0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_A*theta_s + mu_L*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_A/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {
                       // printf("\ninside liquid d2udy2 i = 0 j =1");
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
                else if  (phi_current_uy >  0){
                    // in air
                    if ( phi_next_uy <= 0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_L*theta_n + mu_A*(1-theta_n);
                        d2udy2 = (mu_L/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy <= 0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_L*theta_s + mu_A*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_L/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
                //printf("\nd2udy2[0][1]=%15.10f",d2udy2);
                if (phi_current_ux <= 0) {
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_L + 3.0/Re_L - del*du2dx - del*duvdy );
                }
                else if (phi_current_ux > 0 ){
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_A + 3.0/Re_L - del*du2dx - del*duvdy);
                }
            }
            
            else if (i==nx){
                
                
                u1=(u[i][j]+u[1][j])/2;
                u2=(u[i-1][j]+u[i][j])/2;
                du2dx=(u1*u1-u2*u2)/dx;
                u1=fabs(u[i][j]+u[1][j])*(u[i][j]-u[1][j])/4;
                u2=fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j])/4;
                du2dx += gamma*(u1-u2)/dx;
                //printf("\ndu2dx[nx][1]=%15.10f",du2dx);
                v1=(v[i][j]+v[i+1][j])/2;
                v2=(v[i][j-1]+v[i+1][j-1])/2;
                u1=(u[i][j]+u[i][j+1])/2;
                u2=(u[i][j-1]+u[i][j])/2;
                u3=(u[i][j]-u[i][j+1])/2;
                u4=(u[i][j-1]-u[i][j])/2;
                duvdy=(v1*u1-v2*u2)/dy+gamma*(fabs(v1)*u3-fabs(v2)*u4)/dy;
                //printf("\nduvdy[nx][1]=%15.10f",duvdy);
                //gx =  3/Re_L;
                
                //determine phi for the
                phi_current_ux = (phi[i][j] + phi[i+1][j])/2.0;
                //phi_next_ux    = (phi[i+1][j] + phi[i+2][j])/2.0;
                phi_prev_ux    = (phi[i-1][j] + phi[i][j])/2.0;
                
                phi_current_uy = (phi[i][j] + phi[i+1][j])/2.0;
                phi_next_uy    = (phi[i][j+1] + phi[i+1][j+1])/2.0;
                phi_prev_uy    = (phi[i][j-1] + phi[i+1][j-1])/2.0;
                
                
                // determine d2udx2
                if (phi_current_ux  <=  0){
                    // in liquid
                    /* if ( phi_next_ux>0){
                     theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                     mu_hat = mu_A*theta_e + mu_L*(1-theta_e);
                     d2udx2 = (mu_A/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[i-1][j])/(dx*dx);
                     }*/
                    if  (phi_prev_ux>0) {
                        theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                        mu_hat = mu_A*theta_w + mu_L*(1-theta_w);
                        d2udx2 = (u[1][j]-u[i][j])/(dx*dx)-(mu_A/mu_hat)*(u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else {
                        d2udx2=(u[i-1][j]-2*u[i][j]+u[1][j])/(dx*dx);
                    }
                    
                }
                else if  (phi_current_ux >  0){
                    // in air
                    /*if ( phi_next_ux <= 0){
                     theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                     mu_hat = mu_L*theta_e + mu_A*(1-theta_e);
                     d2udx2 = (mu_L/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[i-1][j])/(dx*dx);
                     }*/
                    if  (phi_prev_ux <= 0){
                        theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                        mu_hat = mu_L*theta_w + mu_A*(1-theta_w);
                        d2udx2 = (u[1][j]-u[i][j])/(dx*dx)-(mu_L/mu_hat)*(u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else {
                        d2udx2=(u[i-1][j]-2*u[i][j]+u[1][j])/(dx*dx);
                    }
                    
                }
                //printf("\n\nd2udx2[nx][1]=%15.10f",d2udx2);
                // determine d2vdy2
                if (phi_current_uy  <=  0){
                    // in liquid
                    if ( phi_next_uy>0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_A*theta_n + mu_L*(1-theta_n);
                        d2udy2 = (mu_A/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy>0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_A*theta_s + mu_L*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_A/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
                else if  (phi_current_uy >  0){
                    // in air
                    if ( phi_next_uy <= 0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_L*theta_n + mu_A*(1-theta_n);
                        d2udy2 = (mu_L/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy <= 0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_L*theta_s + mu_A*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_L/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
                //printf("\n\nd2udy2[nx][1]=%15.10f",d2udy2);
                if (phi_current_ux <= 0) {
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_L + 3.0/Re_L - del*du2dx - del*duvdy );
                }
                else if (phi_current_ux > 0 ){
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_A + 3.0/Re_L - del*du2dx - del*duvdy);
                }
            }
          else  {
                
                u1=(u[i][j]+u[i+1][j])/2;
                u2=(u[i-1][j]+u[i][j])/2;
                du2dx=(u1*u1-u2*u2)/dx;
                u1=fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])/4;
                u2=fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j])/4;
                du2dx += gamma*(u1-u2)/dx;
                //printf("\ndu2dx[i][j]=%15.10f\ti=%d,  j=%d",du2dx,i,j);
                v1=(v[i][j]+v[i+1][j])/2;
                v2=(v[i][j-1]+v[i+1][j-1])/2;
                u1=(u[i][j]+u[i][j+1])/2;
                u2=(u[i][j-1]+u[i][j])/2;
                u3=(u[i][j]-u[i][j+1])/2;
                u4=(u[i][j-1]-u[i][j])/2;
                duvdy=(v1*u1-v2*u2)/dy+gamma*(fabs(v1)*u3-fabs(v2)*u4)/dy;
                //printf("\nduvdy[0][1]=%15.10fti=%d,  j=%d",duvdy,i,j);
                //gx =  3/Re_L;
                //getchar();
                //determine phi for the
                phi_current_ux = (phi[i][j] + phi[i+1][j])/2.0;
                phi_next_ux    = (phi[i+1][j] + phi[i+2][j])/2.0;
                phi_prev_ux    = (phi[i-1][j] + phi[i][j])/2.0;
                
                phi_current_uy = (phi[i][j] + phi[i+1][j])/2.0;
                phi_next_uy    = (phi[i][j+1] + phi[i+1][j+1])/2.0;
                phi_prev_uy    = (phi[i][j-1] + phi[i+1][j-1])/2.0;
                
                
                // determine d2udx2
                if (phi_current_ux  <=  0){
                    // in liquid
                    if ( phi_next_ux>0){
                        theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                        mu_hat = mu_A*theta_e + mu_L*(1-theta_e);
                        d2udx2 = (mu_A/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else if  (phi_prev_ux>0){
                        theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                        mu_hat = mu_A*theta_w + mu_L*(1-theta_w);
                        d2udx2 = (u[i+1][j]-u[i][j])/(dx*dx)-(mu_A/mu_hat)*(u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else {
                        d2udx2=(u[i-1][j]-2*u[i][j]+u[i+1][j])/(dx*dx);
                    }
                    
                }
                else if  (phi_current_ux >  0){
                    // in air
                    if ( phi_next_ux <= 0){
                        theta_e = fabs(phi_current_ux/(phi_current_ux - phi_next_ux));
                        mu_hat = mu_L*theta_e + mu_A*(1-theta_e);
                        d2udx2 = (mu_L/mu_hat)*(u[i+1][j]-u[i][j])/(dx*dx) - (u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else if  (phi_prev_ux <= 0){
                        theta_w = fabs(phi_current_ux/(phi_current_ux - phi_prev_ux));
                        mu_hat = mu_L*theta_w + mu_A*(1-theta_w);
                        d2udx2 = (u[i+1][j]-u[i][j])/(dx*dx)-(mu_L/mu_hat)*(u[i][j]-u[i-1][j])/(dx*dx);
                    }
                    else {
                        d2udx2=(u[i-1][j]-2*u[i][j]+u[i+1][j])/(dx*dx);
                    }
                    
                }
               //printf("\nd2udx2[i][j]=%15.10f\ti=%d,  j=%d",d2udx2,i,j);
                // determine d2vdy2
                if (phi_current_uy  <=  0){
                    // in liquid
                    if ( phi_next_uy>0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_A*theta_n + mu_L*(1-theta_n);
                        d2udy2 = (mu_A/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy>0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_A*theta_s + mu_L*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_A/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
                else if  (phi_current_uy >  0){
                    // in air
                    if ( phi_next_uy <= 0){
                        theta_n = fabs(phi_current_uy/(phi_current_uy - phi_next_uy));
                        mu_hat = mu_L*theta_n + mu_A*(1-theta_n);
                        d2udy2 = (mu_L/mu_hat)*(u[i][j+1]-u[i][j])/(dy*dy) - (u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else if  (phi_prev_uy <= 0){
                        theta_s = fabs(phi_current_uy/(phi_current_uy - phi_prev_uy));
                        mu_hat = mu_L*theta_s + mu_A*(1-theta_s);
                        d2udy2 = (u[i][j+1]-u[i][j])/(dy*dy)-(mu_L/mu_hat)*(u[i][j]-u[i][j-1])/(dy*dy);
                    }
                    else {//printf("\ninside liquid d2udy2 i = %d j =%d",i,j);
                        d2udy2=(u[i][j-1]-2*u[i][j]+u[i][j+1])/(dy*dy);
                    }
                    
                }
               // printf("\nd2udy2[i][j]=%15.10f\ti=%d,  j=%d, phi_current_uy=%15.10f, phi_next_uy=%15.10f,phi_prev_uy=%15.10f",d2udy2,i,j,phi_current_uy,phi_next_uy,phi_prev_uy);
                if (phi_current_ux <= 0) {
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_L + 3.0/Re_L - del*du2dx - del*duvdy );
                }
                else if (phi_current_ux > 0 ){
                    
                    F[i][j] = del*u[i][j] + dt*( (del*del*d2udx2 + d2udy2)/Re_A + 3.0/Re_L - del*du2dx - del*duvdy);
                }
            }
       
           }
         }
    /* Compute G */
    for (i=1; i<=nx; i++) {
        
        /* j=0; j=ny */
        G[i][0]=del*del*v[i][0];
        G[i][ny]=del*del*v[i][ny];
        //G[i][0]= ( dt* (v[i][2]-2.0*v[i][1])/(dy*dy) )*(del/Re_L) -  ((3.0*cos(beta)/sin(beta))/Re_L)*dt;
        //G[i][ny]= del*del*v[i][ny] - ((3.0*cos(beta)/sin(beta))/Re_L)*dt;
        //G[i][ny] = del*del*v[i][ny] + dt*( -(3.0*cos(beta)/sin(beta))/Re_L + del*del*del*(v[i-1][ny] -2.0*v[i][ny] + v[i+1][ny])/(dx*dx*Re_L)
	//	+ del*(v[i][ny-2] - 2.0*v[i][ny-1])/(dy*dy*Re_L) - del*del*(v[i][j]*(u[i][ny]-u[i-1][ny])/dx + u[i][j]*(v[i][ny]-v[i-1][ny])/dx ));
        for (j=1; j<ny; j++) {
            
            v1=(v[i][j]+v[i][j+1])/2;
            v2=(v[i][j-1]+v[i][j])/2;
            dv2dy=(v1*v1-v2*v2)/dy;
            v1=fabs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])/4;
            v2=fabs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j])/4;
            dv2dy += gamma*(v1-v2)/dy;
            
            u1=(u[i][j]+u[i][j+1])/2;
            u2=(u[i-1][j]+u[i-1][j+1])/2;
            v1=(v[i][j]+v[i+1][j])/2;
            v2=(v[i-1][j]+v[i][j])/2;
            v3=(v[i][j]-v[i+1][j])/2;
            v4=(v[i-1][j]-v[i][j])/2;
            duvdx=(u1*v1-u2*v2)/dx+gamma*(fabs(u1)*v3-fabs(u2)*v4)/dx;
            
            gy=0.0;;
            
            //determine phi
            phi_current_vx = (phi[i][j] + phi[i][j+1])/2.0;
            phi_next_vx    = (phi[i+1][j] + phi[i+1][j+1])/2.0;
            phi_prev_vx    = (phi[i-1][j] + phi[i-1][j+1])/2.0;
            
            phi_current_vy = (phi[i][j] + phi[i][j+1])/2.0;
            phi_next_vy    = (phi[i][j+1] + phi[i][j+2])/2.0;
            phi_prev_vy    = (phi[i][j-1] + phi[i][j])/2.0;
            
            
            // determine d2vdx2
            if (phi_current_vx  <=  0){
                // in liquid
                if ( phi_next_vx>0){
                    theta_e = fabs(phi_current_vx/(phi_current_vx - phi_next_vx));
                    mu_hat = mu_A*theta_e + mu_L*(1-theta_e);
                    d2vdx2 = (mu_A/mu_hat)*(v[i+1][j]-v[i][j])/(dx*dx) - (v[i][j]-v[i-1][j])/(dx*dx);
                }
                else if  (phi_prev_vx>0){
                    theta_w = fabs(phi_current_vx/(phi_current_vx - phi_prev_vx));
                    mu_hat = mu_A*theta_w + mu_L*(1-theta_w);
                    d2vdx2 = (v[i+1][j]-v[i][j])/(dx*dx)-(mu_A/mu_hat)*(v[i][j]-v[i-1][j])/(dx*dx);
                }
                else {
                    d2vdx2=(v[i-1][j]-2*v[i][j]+v[i+1][j])/(dx*dx);
                }
                
            }
            else if  (phi_current_vx >  0){
                // in air
                if ( phi_next_vx <= 0){
                    theta_e = fabs(phi_current_vx/(phi_current_vx - phi_next_vx));
                    mu_hat = mu_L*theta_e + mu_A*(1-theta_e);
                    d2vdx2 = (mu_L/mu_hat)*(v[i+1][j]-v[i][j])/(dx*dx) - (v[i][j]-v[i-1][j])/(dx*dx);
                }
                else if  (phi_prev_vx <= 0){
                    theta_w = fabs(phi_current_vx/(phi_current_vx - phi_prev_vx));
                    mu_hat = mu_L*theta_w + mu_A*(1-theta_w);
                    d2vdx2 = (v[i+1][j]-v[i][j])/(dx*dx)-(mu_L/mu_hat)*(v[i][j]-v[i-1][j])/(dx*dx);
                }
                else {
                    d2vdx2=(v[i-1][j]-2*v[i][j]+v[i+1][j])/(dx*dx);
                }
                
            }
            
            // determine d2vdy2
            if (phi_current_vy  <=  0){
                // in liquid
                if ( phi_next_vy>0){
                    theta_n = fabs(phi_current_vy/(phi_current_vy - phi_next_vy));
                    mu_hat = mu_A*theta_n + mu_L*(1-theta_n);
                    d2vdy2 = (mu_A/mu_hat)*(v[i][j+1]-v[i][j])/(dy*dy) - (v[i][j]-v[i][j-1])/(dy*dy);
                }
                else if  (phi_prev_vy>0){
                    theta_s = fabs(phi_current_vy/(phi_current_vy - phi_prev_vy));
                    mu_hat = mu_A*theta_s + mu_L*(1-theta_s);
                    d2vdy2 = (v[i][j+1]-v[i][j])/(dy*dy)-(mu_A/mu_hat)*(v[i][j]-v[i][j-1])/(dy*dy);
                }
                else {
                    d2vdy2=(v[i][j-1]-2*v[i][j]+v[i][j+1])/(dy*dy);
                }
                
            }
            else if  (phi_current_vy >  0){
                // in air
                if ( phi_next_vy <= 0){
                    theta_n = fabs(phi_current_vy/(phi_current_vy - phi_next_vy));
                    mu_hat = mu_L*theta_n + mu_A*(1-theta_n);
                    d2vdy2 = (mu_L/mu_hat)*(v[i][j+1]-v[i][j])/(dy*dy) - (v[i][j]-v[i][j-1])/(dy*dy);
                }
                else if  (phi_prev_vy <= 0){
                    theta_s = fabs(phi_current_vy/(phi_current_vy - phi_prev_vy));
                    mu_hat = mu_L*theta_s + mu_A*(1-theta_s);
                    d2vdy2 = (v[i][j+1]-v[i][j])/(dy*dy)-(mu_L/mu_hat)*(v[i][j]-v[i][j-1])/(dy*dy);
                }
                else {
                    d2vdy2=(v[i][j-1]-2*v[i][j]+v[i][j+1])/(dy*dy);
                }
                
            }
	    // printf("\ngravity term = %15.10f", (3*cos(beta)/sin(beta)));
            if (phi_current_vx <= 0) {
                
                G[i][j] = del*del*v[i][j] + dt*( (del*del*del*d2vdx2 + del*d2vdy2)/Re_L - (3*cos(beta)/sin(beta))/Re_L - del*del*duvdx - del*del*dv2dy );
            }
            else if (phi_current_vx > 0 ){
                
                G[i][j] = del*del*v[i][j] + dt*( (del*del*del*d2vdx2 + del*d2vdy2)/Re_A - (3*cos(beta)/sin(beta))/Re_L - del*del*duvdx - del*del*dv2dy );
            }
            
        }
    }
    
    /* Compute H */
    for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            H[i][j]=(del*(F[i][j]-F[i-1][j])/dx+(G[i][j]-G[i][j-1])/dy)/dt;
        }
    }
    
    return 0;
}
