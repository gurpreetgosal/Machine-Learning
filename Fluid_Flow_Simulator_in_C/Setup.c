#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>

/*
 Read parameters nx, ny, time from file parm.
 Allocate space for u, v, p.
 */

int Setup(char *infile,int *nx,int *ny,double *time,double ***u,double ***v,double ***p)
{
    int    i,Nu,Nv,Np,Nphi,num;
    double **U,**V,**P;
    char   buf[20];
    FILE   *f;
    printf(" inside Setup\n");
    
    /* Read from file */
    if ((f=fopen(infile,"r"))==NULL) {
        printf("Cannot open file %s\n",infile);
        exit(0);
    }
    fscanf(f,"nx=%d\n",nx);
    fscanf(f,"ny=%d\n",ny);
    fscanf(f,"time=%s\n",buf);
    *time = atof(buf);
    fclose(f);
    
    
    /* Allocate space for u, v, p */
    Nu   = (*nx+1)*(*ny+2);
    Nv   = (*nx+2)*(*ny+1);
    Np   = (*nx+2)*(*ny+2);
    // Nphi = (*nx+2)*(*ny+2);

    U = (double **) malloc((*nx+1)*sizeof(double *));
    U[0] = (double *) malloc(Nu*sizeof(double));
    for (i=1; i<*nx+1; i++) U[i]=U[i-1]+(*ny+2);
    for (i=0; i<Nu; i++) U[0][i]=0.0;
    
    V = (double **) malloc((*nx+2)*sizeof(double *));
    V[0] = (double *) malloc(Nv*sizeof(double));
    for (i=1; i<*nx+2; i++) V[i]=V[i-1]+(*ny+1);
    for (i=0; i<Nv; i++) V[0][i]=0.0;
    
    P = (double **) malloc((*nx+2)*sizeof(double *));
    P[0] = (double *) malloc(Np*sizeof(double));
    for (i=1; i<*nx+2; i++) P[i]=P[i-1]+(*ny+2);
    for (i=0; i<Np; i++) P[0][i]=0.0;
      
    
    *u=U; *v=V; *p=P;
    printf(" inside Setup finished\n");

    return 0;
}

