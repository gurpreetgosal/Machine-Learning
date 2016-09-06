#include <stdlib.h>
#include <math.h>
#include <stdio.h>


double Norm(int ,int ,double **);
int Residual(int ,int,double **,double **,double **,int,double**,double, double,double);
int matvecprodCG(int,int,double **,double **,double **,double,double,double);
int copy_dir_r(int,int,double **,double **,double **);
double compute_a(int,int ,double **,double **,double **);
int update_p(int, int, double **, double **, double **,double);
int update_residual(int, int, double **, double **, double **,double);
double ratio_r_newr(int, int, double **, double **);
int update_dir(int, int, double **, double **,double);
int copy_r(int, int, double **, double **);
double product_r(int, int, double **);
int set_BC(int, int, double ** );
int p_BC(int, int, double ** );
int p_scale(int nx,int ny,double **, double, double, double **);
double H_norm(int, int , double **);
/*
Solve Pressure equation using CG method
*/
int Pressure(int nx,int ny,double **H,double **p,int counter,double **phi,double rho_L,double rho_A,double del)
{ printf("\ninside pressure\n");
  int    i,j,k=0,maxiter=20000,Np;
  double tol=1e-5,**r,**new_r,**dir,**vprod,a_ratio, r_ratio, rprod, rnorm;
     
  // CG iterations 
  //x = initial guess for inv(A)*b
  //r = b - A*x   ... residual, to be made small
  Np = (nx+2)*(ny+2);
  
  r = (double **) malloc((nx+2)*sizeof(double *));
  r[0] = (double *) malloc(Np*sizeof(double));
  for (i=1; i<nx+2; i++) r[i]=r[i-1]+(ny+2);
  for (i=0; i<Np; i++) r[0][i]=0.0;
  // printf("\nprinting initialized r = %15.10f\n",r[1][1]);
  new_r = (double **) malloc((nx+2)*sizeof(double *));
  new_r[0] = (double *) malloc(Np*sizeof(double));
  for (i=1; i<nx+2; i++) new_r[i]=new_r[i-1]+(ny+2);
  for (i=0; i<Np; i++) new_r[0][i]=0.0;
    
  dir = (double **) malloc((nx+2)*sizeof(double *));
  dir[0] = (double *) malloc(Np*sizeof(double));
  for (i=1; i<nx+2; i++) dir[i]=dir[i-1]+(ny+2);
  for (i=0; i<Np; i++) dir[0][i]=0.0;
   printf("\nprinting initialized dir = %15.10f\n",dir[1][1]);
  vprod = (double **) malloc((nx+2)*sizeof(double *));
  vprod[0] = (double *) malloc(Np*sizeof(double));
  for (i=1; i<nx+2; i++) vprod[i]=vprod[i-1]+(ny+2);
  for (i=0; i<Np; i++) vprod[0][i]=0.0;
 
  p_BC(nx,ny,p);
  Residual(nx,ny,H,p,r,counter,phi,rho_L,rho_A,del);
  printf("\nprinting  r after first calculation residual = %15.10f\n",r[1][1]);
  // p = r         ... initial "search direction"
  copy_dir_r(nx,ny,r,dir,p); // check for validity of this assignment
  set_BC(nx,ny,dir);
  //set update ghost values Neumann BCs
 
  rnorm = product_r(nx,ny,r);
  //until ( new_r'*new_r small enough )
    
  printf("printing phi for current iteration");
  /* for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            printf("\t%15.10f\n",phi[i][j]);
       }
       }*/

  while ( (rnorm>=tol*H_norm(nx,ny,H))&&(k<maxiter))  {
    //printf("\ninside loop\n\n\n");
    k= k+1;
    //getchar();
     
    //printf("\tcurrent time pressure iteration number = %d\n",k);
    //printf("\nrnorm = %15.10f\n",rnorm);
   
    // set the BCs on dir and p
    // set_BC(nx,ny,dir);

    //v = A*p... matrix-vector multiply 
    matvecprodCG(nx,ny,dir,vprod,phi,rho_L,rho_A,del);
    
    //a = ( r'*r ) / ( p'*v )... dot product (BLAS1)
    a_ratio = compute_a(nx,ny,r,dir,vprod);
    //printf("\n \t a_ratio = %15.10f\n\n  ",a_ratio);
   
    //for (i=1;i<=nx;i++){
     
    //printf("\n ri = %15.10f\t  pi = %15.10f\t dir_i = %15.10f\t vprod_i = %15.10f\tHi =%15.10f\t\n",r[1][1],p[1][1],dir[1][1],vprod[1][1],H[1][1]);
       // }
    // printf("\na_ratio = %15.10f\n",a_ratio);
    //x = x + a*p... compute updated solution
    update_p(nx,ny,p,dir,r,a_ratio);
    p_BC(nx,ny,p);
    //new_r = new_r - a*v... compute updated residual
    update_residual(nx,ny,r,new_r,vprod,a_ratio);
    //printf("\nupdated residual\n") ;
     
    //g = ( new_r'*new_r ) / ( r'*r )... dot product (BLAS1)
    r_ratio = ratio_r_newr(nx,ny,r,new_r);
    //printf("\nr_ratio = %15.10f\n",r_ratio);
        
    //p = new_r + g*p... compute updated search direction
    update_dir(nx,ny,new_r,dir,r_ratio);
    set_BC(nx,ny,dir);
    //printf("\ndir updated\n");    
    //r = new_r
    copy_r(nx,ny,r,new_r);
    //printf("\n newr_i = %15.10f\t updated ri =%15.10f \n",new_r[1][1],r[1][1]);  
    // computer r'*r for convergence
    rnorm = product_r(nx,ny,r);
    // printf("\n rnorm = %15.10f",rnorm);
    // break; 
   
  }
  set_BC(nx,ny,dir);
  printf("\npressure iterations= %d\tstopping criteria=%15.10f\trnorm=%15.10f\n",k,tol*H_norm(nx,ny,H),rnorm);
  // printf("pressure for current iteration\n\n");
  //p_scale(nx,ny,p,rho_L,rho_A,phi);
  /*for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            printf("%15.10f,",p[i][j]);
           }
	   }*/

  // getchar();
  //printf("/n iterations = %15.10f",k);
  free(r[0]);  free(r);  free(new_r[0]);  free(new_r);
  free(dir[0]);  free(dir);  free(vprod[0]);  free(vprod); 

  return 0;
}

double Norm(int nx,int ny,double **x)
{
  int    i,j;
  double sum=0.0,nrm;
    
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      sum += x[i][j]*x[i][j];
    }
  }
  nrm = (sum);
    
  return nrm;
}

int Residual(int nx,int ny,double **H,double **p,double **r,int counter,double **phi,double rho_L,double rho_A,double del)
{
   int    i,j;
    double hx,hy;
  //double phi_e,phi_w,phi_n,phi_s,phi_c;
    
  hx=5.0/nx; hy=2.0/ny;
 
    
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      
        double rho_in=0.0, rho_out=0.0,theta_e=0, theta_w=0,theta_s=0,theta_n=0,rho_ave=0,AE=0,AW=0,AS=0,AN=0,AC=0;
            
            if (phi[i][j] <=  0){
	      //printf("\n phi[i][j]<=0 ");
                rho_in  = rho_L;
                rho_out = rho_A;
            }
	    else if (phi[i][j] >  0) {
	      //printf("\n phi[i][j] > 0 ");
                rho_in  = rho_A;
                rho_out = rho_L;
	    }
            //printf("\nphi < 0 inside liquid\n ");
            //AE=(1.0/(hx*hx))*(1.0/rho_in); AW=AE; AS=1.0/(hy*hy)*(1.0/rho_in); AN=AS; AC=-(AE+AW+AN+AS);
            
            //Interface between x_k and x_k+1
            if (phi[i+1][j]>0 != phi[i][j]>0){
                theta_e = fabs(phi[i][j]/(phi[i][j] - phi[i+1][j]));
                rho_ave = rho_in*theta_e + rho_out*(1-theta_e);
                AE = del*del*(1.0/(hx*hx))*(1.0/rho_ave);}
	    else {
               AE=del*del*(1.0/(hx*hx))*(1.0/rho_in);
            }
            //Interface between x_k and x_k-1
            if (phi[i-1][j]>0 != phi[i][j]>0){
                theta_w = fabs(phi[i][j]/(phi[i][j] - phi[i-1][j]));
                rho_ave = rho_in*theta_w + rho_out*(1-theta_w);
                AW =del*del*(1.0/(hx*hx))*(1.0/rho_ave);}
             else {
               AW=del*del*(1.0/(hx*hx))*(1.0/rho_in);
            }
            
            
            //Interface between y_k and y_k+1
           if (phi[i][j+1]>0!= phi[i][j]>0){
                theta_n = fabs( phi[i][j]/(phi[i][j] - phi[i][j+1]));
                rho_ave = rho_in*theta_n + rho_out*(1-theta_n);
                AN = (1.0/(hy*hy))*(1.0/rho_ave); }
            else {
                AN=(1.0/(hy*hy))*(1.0/rho_in);
            
            }
            //Interface between y_k and y_k-1
            if (phi[i][j-1]>0!= phi[i][j]>0){
                theta_s = fabs(phi[i][j]/(phi[i][j] - phi[i][j-1]));
                rho_ave = rho_in*theta_s + rho_out*(1-theta_s);
                AS =(1.0/(hy*hy))*(1.0/rho_ave); }
            else {
                AS=(1.0/(hy*hy))*(1.0/rho_in);   
            }
            
            AC = -(AE + AW + AS + AN);
           

       //printf("\ncheck");
       //printf("\n  %15.10f\t  %15.10f\t %15.10f\t %15.10f\n",theta_e,theta_w,theta_n,theta_s);
       //printf("\n  %e\t  %15.10f\t %15.10f\t AE = %15.10f\tAW =%15.10f \n",phi[i][j],phi[i+1][j],phi[i-1][j],AE,AW);
                 
        r[i][j]=H[i][j]-(AE*p[i+1][j]+AW*p[i-1][j]+AN*p[i][j+1]+AS*p[i][j-1]+AC*p[i][j]);
        
    }
  }
  // getchar();
  return 0;
}


int matvecprodCG(int nx,int ny,double **dir,double **vprod,double **phi, double rho_L,double rho_A,double del){
  int    i,j;
  double hx,hy;
  //double phi_e,phi_w,phi_n,phi_s,phi_c;
    
  hx=5.0/nx; hy=2.0/ny;
 
    for (j=1; j<=ny; j++) {
        for (i=1; i<=nx; i++) {
            
            double rho_in=0.0, rho_out=0.0,theta_e=0, theta_w=0,theta_s=0,theta_n=0,rho_ave=0,AE=0,AW=0,AS=0,AN=0,AC=0;
            
            if (phi[i][j] <=  0){
	      //printf("\n phi[i][j]<=0 ");
                rho_in  = rho_L;
                rho_out = rho_A;
            }
	    else if (phi[i][j] >  0) {
	      //printf("\n phi[i][j] > 0 ");
                rho_in  = rho_A;
                rho_out = rho_L;
	    }
            //printf("\nphi < 0 inside liquid\n ");
            //AE=(1.0/(hx*hx))*(1.0/rho_in); AW=AE; AS=1.0/(hy*hy)*(1.0/rho_in); AN=AS; AC=-(AE+AW+AN+AS);
            
            //Interface between x_k and x_k+1
            if (phi[i+1][j]>0 != phi[i][j]>0){
                theta_e = fabs(phi[i][j]/(phi[i][j] - phi[i+1][j]));
                rho_ave = rho_in*theta_e + rho_out*(1-theta_e);
                AE = del*del*(1.0/(hx*hx))*(1.0/rho_ave);}
	    else {
               AE=del*del*(1.0/(hx*hx))*(1.0/rho_in);
            }
            //Interface between x_k and x_k-1
            if (phi[i-1][j]>0 != phi[i][j]>0){
                theta_w = fabs(phi[i][j]/(phi[i][j] - phi[i-1][j]));
                rho_ave = rho_in*theta_w + rho_out*(1-theta_w);
                AW =del*del*(1.0/(hx*hx))*(1.0/rho_ave);}
             else {
               AW=del*del*(1.0/(hx*hx))*(1.0/rho_in);
            }
            
            
            //Interface between y_k and y_k+1
           if (phi[i][j+1]>0!= phi[i][j]>0){
                theta_n = fabs( phi[i][j]/(phi[i][j] - phi[i][j+1]));
                rho_ave = rho_in*theta_n + rho_out*(1-theta_n);
                AN = (1.0/(hy*hy))*(1.0/rho_ave); }
            else {
                AN=(1.0/(hy*hy))*(1.0/rho_in);
            
            }
            //Interface between y_k and y_k-1
            if (phi[i][j-1]>0!= phi[i][j]>0){
                theta_s = fabs(phi[i][j]/(phi[i][j] - phi[i][j-1]));
                rho_ave = rho_in*theta_s + rho_out*(1-theta_s);
                AS =(1.0/(hy*hy))*(1.0/rho_ave); }
            else {
                AS=(1.0/(hy*hy))*(1.0/rho_in);   
            }
            
            AC = -(AE + AW + AS + AN);
            /*
	      printf("\n\n checking modified matrix terms \n");
	      printf("\n hx= %15.10f\t  hy=%15.10f\t  rho_in=%15.10f\t rho_ave=%15.10f\t",hx,hy,rho_in,rho_ave);
	    printf("\n  %e\t  %15.10f\t %15.10f\t AE = %15.10f\tAW =%15.10f\t AC =%15.10f ",phi[i][j],phi[i+1][j],phi[i-1][j],AE,AW,AC);
	    printf("\n  %e\t  %15.10f\t %15.10f\t AN = %15.10f\tAS =%15.10f\t AC =%15.10f \n",phi[i][j],phi[i][j+1],phi[i][j-1],AN,AS,AC);
	    */
        vprod[i][j]= (AE*dir[i+1][j]+AW*dir[i-1][j]+AN*dir[i][j+1]+AS*dir[i][j-1]+AC*dir[i][j]);
    }
  }
    return 0;
}


int copy_dir_r(int nx,int ny,double **r,double **dir,double **p){
  int i,j;
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      dir[i][j] = r[i][j];
     
    }
  }
  return 0;
}

double compute_a(int nx,int ny,double **r,double **dir,double **vprod){
  int i,j;
  double sum_rr=0.0,sum_dirvprod =0.0, a_ratio;
    
  //a = ( r'*r ) / ( p'*v )
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      sum_rr = sum_rr + r[i][j]*r[i][j];
      sum_dirvprod = sum_dirvprod + dir[i][j]*vprod[i][j];
    }
  }
  a_ratio = sum_rr/sum_dirvprod;
  
  return a_ratio;
}

int update_p(int nx, int ny, double **p, double **dir, double **r,double a_ratio){
  int i,j;
  //x = x + a*p... compute updated solution
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      p[i][j] += a_ratio*dir[i][j];
    }
  }
  // printf("\nwithin the p update function\n");
  //printf("\n%15.10f\t%15.10f\t%15.10f\n",r[i][j],dir[i][j],p[i][j]);
  return 0;
}

int update_residual(int nx, int ny, double **r, double **new_r, double **vprod, double a_ratio){
  int i,j;
  //new_r = r - a*v... compute updated residual
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      new_r[i][j] = r[i][j] - a_ratio*vprod[i][j];
    }
  }
  return 0;
}


double ratio_r_newr(int nx, int ny, double **r, double **new_r ){
  //g = ( new_r'*new_r ) / ( r'*r )... dot product (BLAS1)
  int i,j;
  double sum_r=0.0, sum_newr=0.0,r_ratio;
    
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      sum_r = sum_r + r[i][j]*r[i][j];
      sum_newr = sum_newr + new_r[i][j]*new_r[i][j];        
    }
  }
  r_ratio = sum_newr/sum_r; 
  return r_ratio;
    
}

int update_dir(int nx, int ny, double **new_r, double **dir, double r_ratio){
  int i,j;
  //p = new_r + g*p
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      dir[i][j] = new_r[i][j]+ r_ratio*dir[i][j];
    }
  }
  return 0;
}

int copy_r(int nx, int ny, double **r, double **new_r){
  int i,j;
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      r[i][j] = new_r[i][j];
    }
  }
  return 0;
}

double product_r(int nx, int ny, double **r){
  int i,j;
  double sum_rprod=0.0,rnorm,rprod;
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      sum_rprod = sum_rprod + r[i][j]*r[i][j];
    }
  }
  rnorm = sum_rprod;
  rnorm = sqrt(sum_rprod);
  return rnorm;
}

int set_BC(int nx, int ny, double **dir ){
  int i,j;
  /* set ghost values*/
  /* for (i=0; i<=nx+1; i++) {
    dir[i][0]=dir[i][1];
    dir[i][ny+1]=dir[i][ny];
  }
  for (j=0; j<=ny+1; j++) {
    dir[0][j]=dir[1][j];
    dir[nx+1][j]=dir[nx][j];
    }*/
  for (i=0; i<=nx+1; i++) {
    dir[i][0]=dir[i][1];
    dir[i][ny+1]=dir[i][ny];
    }
  for (j=0; j<=ny+1; j++) {
     dir[0][j]=dir[nx-1][j];
     dir[1][j]=dir[nx][j];
     dir[nx+1][j]=dir[2][j];
     } 
  //dir[0][0]=dir[1][1];
  //dir[nx+1][0]=dir[nx][1];
  //dir[0][ny+1]=dir[1][ny];
  //dir[nx+1][ny+1]=dir[nx][ny];
  
  return 0;
}

int p_BC(int nx, int ny, double **p ){
  int i,j;
  /* set ghost values*/
  /*  for (i=0; i<=nx+1; i++) {
    p[i][0]=p[i][1];
    p[i][ny+1]=p[i][ny];
  }
  for (j=0; j<=ny+1; j++) {
    p[0][j]=p[1][j];
    p[nx+1][j]=p[nx][j];
    }*/
  for (i=0; i<=nx+1; i++) {
    p[i][0]=p[i][1];
    p[i][ny+1]=p[i][ny];
  }
  for (j=0; j<=ny+1; j++) {
    p[0][j]=p[nx-1][j];
    p[nx+1][j]=p[2][j];
    p[1][j]=p[nx][j];
    
    }

  //p[0][0]=p[1][1];   p[nx+1][0]=p[nx][1];
  //p[0][ny+1]=p[1][ny]; p[nx+1][ny+1]=p[nx][ny];
    
  return 0;
}

int p_scale(int nx, int ny, double **p, double rho_L, double rho_A, double **phi ){
  int i, j;
  double rho_in;
   
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
     
      if (phi[i][j] <=  0){   
          rho_in  = rho_L;
      }
      else if (phi[i][j] >  0) {
          rho_in  = rho_A;
      }
      /* printf("\n pressure scaling");
      printf("\n before scalin rho_in = %15.10f", rho_in);
      printf("\n after scaling p = %15.10f",p[i][j]);*/
      //printf("\n  rho_in = %15.10f\t p[i][j]= %15.10f\t p[i][j]/rho_in = %15.10f ", rho_in,p[i][j],p[i][j]/rho_in);
      p[i][j]= p[i][j]/rho_in;
      //printf("\n  after after rho_in = %15.10f\t p[i][j]= %15.10f\t p[i][j]/rho_in = %15.10f ", rho_in,p[i][j],p[i][j]/rho_in);
    }
  }
    return 0;
  
}
double H_norm(int nx, int ny, double **H){
  int i,j;
  double sum_H=0.0,Hnorm,Hprod;
  for (j=1; j<=ny; j++) {
    for (i=1; i<=nx; i++) {
      sum_H = sum_H + H[i][j]*H[i][j];
      //printf("\nH=%15.10f",H[i][j]);
    }
  }
  //Hnorm = sum_rprod;
  //printf("\nsum_H=%15.10f",sum_H);
  Hnorm = sqrt(sum_H);
  return Hnorm;
}
