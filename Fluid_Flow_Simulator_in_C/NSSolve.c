#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
/*
 Solve the Navier-Stokes equations and
 output the results to output.m
 */

int NSSolve(int nx,int ny,double time,double **u,double **v,double **p)
{
    
    printf("inside NSSOlve");
    int    k,nt,i,j,counter=0,Nphi;
    double **F,**G,**H,**phi,sum_H,h0,x_current,y_current,x_dim = 5.0,y_dim=2.0;
    double t,dt,x,y,del_x,del_y,xc,yc, rho_L = 1000.0,rho_A = 1.0, t_store[100];
    double del = 0.1, mu_L = 1000.0/18.0, mu_A= 1.0, beta = M_PI/4.0 ;
    //double Re_L = rho_L*1.0*1.0/mu_L , Re_A = rho_A*1.0*1.0/mu_A;
    double Re_L = 1.0,  Re_A = 1.0/15.0;
    int  t_count = 1, t_writecount[100];
    FILE   *f,*g;
    del_x = 5.0/nx;
    del_y = 2.0/ny;
    Nphi = (nx+2)*(ny+2);
    
    /* Initialization */
    nt=0;
    t=0.0;
    dt = 0.02;
    
    F = (double **) malloc((nx+1)*sizeof(double *));
    F[0] = (double *) malloc((nx+1)*(ny+2)*sizeof(double));
    for (k=1; k<nx+1; k++) F[k]=F[k-1]+ny+2;
    G = (double **) malloc((nx+2)*sizeof(double *));
    G[0] = (double *) malloc((nx+2)*(ny+1)*sizeof(double));
    for (k=1; k<nx+2; k++) G[k]=G[k-1]+ny+1;
    H = (double **) malloc((nx+2)*sizeof(double *));
    H[0] = (double *) malloc((nx+2)*(ny+2)*sizeof(double));
    for (k=1; k<nx+2; k++) H[k]=H[k-1]+ny+2;
    //work = (double *) malloc(4*(nx+2)*(ny+2)*sizeof(double));
    
    phi = (double **) malloc((nx+2)*sizeof(double *));
    phi[0] = (double *) malloc(Nphi*sizeof(double));
    for (i=1; i<nx+2; i++) phi[i]=phi[i-1]+(ny+2);
    for (i=0; i<Nphi; i++) phi[0][i]=0.0;
    
    // fprintf(f,"phi=[\n");
    if ((g=fopen("output.m","w"))==NULL) {
        printf("Cannot open file %s\n","output.m");
        exit(0);
    }
    
    // Initialize phi
    
    for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            
            x = (i-1)*(del_x)+ del_x/2;
            y = (j-1)*(del_y)+ del_y/2;
            xc = 0.5;
            yc = 0.5;
            // phi[i][j] = pow((x-xc),2.0) + pow((y-yc),2.0) - pow(del_x*0.5432,2);
	    //phi[i][j] = sqrt(pow( ( (i*1.0/(nx*1.0)) - 0.5 ),2)+ pow( ( (j*1.0/(ny*1.0)) - 0.5 ),2)) -0.157;
	    //phi[i][j] =  ( (j-1)*del_y + del_y/2.0 ) - 1.0 - 0.1*sin(2.0*M_PI*((i-1)*del_x+del_x/2.0)/10.0) ;
            phi[i][j] = ( (j-1)*del_y + del_y/2.0 ) - 1.0;
            // getchar();
            
        }
    }
    
    
    for (i=0; i<=nx+1; i++) {
        phi[i][0]= phi[i][1];
        phi[i][ny+1]=phi[i][ny];
    }
    for (j=0; j<=ny+1; j++) {
        phi[0][j]=phi[nx-1][j];
        phi[nx+1][j]=phi[2][j];
        phi[1][j]=phi[nx][j];
    }
    /*phi[0][0]=phi[nx][ny];
    phi[nx+1][0]=phi[1][ny];
    phi[0][ny+1]=phi[nx][1];
    phi[nx+1][ny+1]=phi[1][1];*/
    // printf("\nphi initialized");
    
    /* for (j=1; j<=ny; j++) {
       for (i=0; i<=nx; i++) {
         if  ((phi[i][j] + phi[i+1][j])/2.0 <=0) {
           y_current = (j-1)*(del_y+del_y/2.0);
	   x_current = i*del_x;
	   //h0 = 1.0 + 0.1*sin(2*M_PI*i*del_x/10.0);
	   //u[i][j] =( 2.0/(3.0*h0*h0*h0)) * (2*h0 - ((j-1)*del_y+del_y/2.0)) * ((j-1)*del_y+del_y/2.0);  
           // u[i][j] = (1.5)*((j-1)*del_y+del_y/2.0)*(2.0 - ((j-1)*del_y+del_y/2.0) ) ;
           u[i][j] = 1.5*(2.0*y_current-y_current*y_current)*( 1.0 + 0.0*sin(2*M_PI*x_current/x_dim) );         
	 }
	 else{
           //u[i][j] = 0.0; 
            u[i][j] = 1.5;  
            }
       }
       }
     for (j=0; j<ny; j++) {
       for (i=1; i<=nx; i++) {
         if  ((phi[i][j] + phi[i][j+1])/2.0 <=0) {
           y_current = j*del_y;
	   x_current = (i-1)*(del_x+del_x/2.0);
	   //h0 = 1.0 + 0.1*sin(2*M_PI*i*del_x/10.0);
	   //u[i][j] =( 2.0/(3.0*h0*h0*h0)) * (2*h0 - ((j-1)*del_y+del_y/2.0)) * ((j-1)*del_y+del_y/2.0);  
          // u[i][j] = (1.5)*((j-1)*del_y+del_y/2.0)*(2.0 - ((j-1)*del_y+del_y/2.0) ) ;
	 v[i][j] = (0.0*M_PI/x_dim)*(3.0*y_current*y_current -y_current*y_current*y_current )* (cos(2*M_PI*x_current/x_dim));         
	 }
	 else{
           //u[i][j] = 0.0; 
            v[i][j] = 0.0;  
            }
       }
       }
    //printf("\n after u init u11=%15.10f",u[1][1]);
    BC(nx,ny,u,v);*/
    
    //printf("\noutside while");
    /* Evolution */
    
    // specify the time points where results are to be stored
    //double val_c = 0.0;
    *t_store = 0.0;
    
    while (*(t_store+t_count-1) < time){
      *(t_store+t_count) += *(t_store + t_count-1) + 0.1;  
     
       printf("t_store_%d = %15.10f\n ",t_count,*(t_store+t_count));
        t_count++;
       }
    //getchar();
    
    while (t<time-1e-14) {
    //printf("\ninside while");
        counter = counter +1;
        // printf("\n------------time counter = %d------------------\n",counter);
        if (t>=i*time/10-1e-12) {
            printf("Current time = %15.10f\n",t);
            i++;
        }
         
        TimeStep(nx,ny,u,v,&dt,Re_L,Re_A,del);
        
        RHS(nx,ny,dt,u,v,F,G,H,phi,del,rho_L,rho_A,mu_L,mu_A,beta);
	// RHS_check(nx,ny,dt,u,v,F,G,H);
        
        //printf("\npressure ");
	/* if (counter%10000==0){
	   getchar();}*/
        sum_H=0.0;
	 //	flag_H =0 ;
        for (k=1; k<=ny; k++) {
            for (i=1; i<=nx; i++) {
	      //sum_H += fabs(H[i][k]);
	      sum_H += fabs(H[i][k]); 
            } 
            //printf("\n");
        }
	/* printf("\nH\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",H[i][k]);
                
            }
            printf("\n");
	    }*/
        //printf("\nsum_H = %15.10f",sum_H);
	if (sum_H==0.0){
            for (k=1; k<=ny; k++) {
               for (i=1; i<=nx; i++) {
                 p[i][k]=0.0;           
                 }
            //printf("\n");
              }
         }
        else{
              Pressure(nx,ny,H,p,counter,phi,rho_L,rho_A,del);
	      }
        //Pressure(nx,ny,H,p,counter,phi,rho_L,rho_A,del);
        //Pressure(nx,ny,H,p,counter,phi,rho_L,rho_A);
        printf("\n dt= %15.10f",dt);
        
	  printf("\nH\n");
	  /* for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",H[i][k]);
                
            }
            printf("\n");
        }
        printf("\n\n");
	 printf("\nF\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",F[i][k]);
                
            }
            printf("\n");
        }
        printf("\n\n");
        printf("\nG\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",G[i][k]);
                
            }
            printf("\n");
        }
        printf("\n\n");
         printf("\npressure\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",p[i][k]);
                
            }
            printf("\n");
        }
        printf("\n\nphi\n\n");*/
        
        
        // print stuff to file
	//	printf("\n integer t = %d",(int) t);
	if ( counter%50000==0 ||  counter==1  ||  counter==20  ||  counter==50 ||  counter==100  ||  counter==200    ) {
	  // *(t_writecount+t_count) = 1 ;
	  //t_count++;
        fprintf(g,"Current_time = %15.10f\n",t);
         printf("\nphi_inside\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<nx+1; i++) {
                printf("%15.10f,",phi[i][k]);
                
            }
            printf("\n");
        }
        //printf("\n\nphi\n\n");
        
        fprintf(g,"u%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",u[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
        
        fprintf(g,"v%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",v[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
      
        /* fprintf(g,"uu%d=[\n",counter);
     for (k=1; k<=ny; k++) {
     for (i=0; i<nx+1; i++) {
     fprintf(g,"%15.10f\n",u[i][k]);
     }
     }
          fprintf(g,"];\n");

     fprintf(g,"vv%d=[\n",counter);
     for (k=0; k<ny+1; k++) {
     for (i=1; i<=nx; i++) {
     fprintf(g,"%15.10f\n",v[i][k]);
     }
     }
     fprintf(g,"];\n");
       
       
        
        fprintf(g,"F%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",F[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
        
        fprintf(g,"G%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",G[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
        
        fprintf(g,"H%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",H[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");*/
        
        fprintf(g,"p%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",p[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");

        fprintf(g,"phi%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",phi[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
        }
        
       
	printf("u velocity\n\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<=nx; i++) {
                printf("%15.10f,",u[i][k]);
                
            }
            printf("\n");
	    }
	/* printf("v  velocity\n\n");
        for (k=1; k<=ny; k++) {
            for (i=1; i<=nx; i++) {
                printf("%15.10f,",v[i][k]);
                
            }
            printf("\n");
	    }*/

        Velocity(counter,nx,ny,dt,F,G,p,u,v,phi,del,rho_L,rho_A,Re_L);
        BC(nx,ny,u,v);
        printf("\n------------time counter = %d------------------\n",counter);
        printf("Current time = %15.10f\n",t);
        //printf("\n p[1][1]=%15.10f",p[1][1]);
        //getchar();
        
        level_set(nx,ny,dt,u,v,phi);
        
        nt++;
        t+=dt;
    }
    
    for (k=1; k<=ny; k++) {
        for (i=0; i<nx+1; i++) {
            printf("\npressure_final = %15.10f\n",p[i][k]);
        }
    }
    
    /*fprintf(f,"nx=%d; ny=%d;\n",nx,ny);
     fprintf(f,"u=[\n");
     for (k=1; k<=ny; k++) {
     for (i=0; i<nx+1; i++) {
     fprintf(f,"%15.10f\n",u[i][k]);
     }
     }
     fprintf(f,"];\n");
     fprintf(f,"v=[\n");
     for (k=0; k<ny+1; k++) {
     for (i=1; i<=nx; i++) {
     fprintf(f,"%15.10f\n",v[i][k]);
     }
     }
     fprintf(f,"];\n");*/
     fprintf(g,"phi%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",phi[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
     fprintf(g,"u%d=[\n",counter);
        for (k=1; k<=ny; k++) {
            for (i=0; i<nx+1; i++) {
                fprintf(g,"%15.10f,",u[i][k]);
            }
            fprintf(g,"\n");
        }
        fprintf(g,"];\n");
    fprintf(g,"\n counter=%d",counter);
    fprintf(g,"\n nx=%d; ny=%d;\n",nx,ny);
    fprintf(g,"Current_time = %15.10f\n",t);
    //fclose(f);
    fclose(g);
    
    
    free(F[0]); free(F); free(G[0]); free(G); free(H[0]); free(H);
    free(phi[0]); free(phi);
    
    return 0;
}

