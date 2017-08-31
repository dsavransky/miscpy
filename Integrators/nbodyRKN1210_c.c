/*=================================================================
 *
 * NBODYRKN1210_C.C	RKN for NBODY Equation
 *
 * NBODYRKN_C(T_IN, Y_IN, DY_IN, MUS, n, Y_OUT, DY_OUT)
 *      T_IN    2  x 1 : [t_0, t_f]
 *      Y_IN    3n x 1 : initial positions
 *      DY_IN   3n x 1 : initial velocities
 *      MUS     n  x 1 : Gravitational paramters
 *	n	scalar : number of bodies
 *
 *=================================================================*/
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "RKN1210consts.h"
#include "nbodyRKN1210_c.h"

#if !defined(NORM)
#define NORM(x1, x2, x3) sqrt(pow(x1, 2)+pow(x2, 2)+pow(x3, 2))
#endif

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

double arrMax(double a[], int size) {
    double maxVal = a[0];
    int i;
    for (i=1; i<size; i++) {
        if (a[i] > maxVal) { maxVal = a[i]; }
    }
    return maxVal;
}

double arrMin(double a[], int size) {
    double minVal = a[0];
    int i;
    for (i=1; i<size; i++) {
        if (a[i] < minVal) { minVal = a[i]; }
    }
    return minVal;
}

double arrAbsMax(double a[], int size) {
    double maxVal = fabs(a[0]);
    int i;
    for (i=1; i<size; i++) {
        if (fabs(a[i]) > maxVal) { maxVal = fabs(a[i]);}
    }
    return maxVal;
}

double arrAbsMin(double a[], int size) {
    double minVal = fabs(a[0]);
    int i;
    for (i=1; i<size; i++) {
        if (fabs(a[i]) < minVal) { minVal = fabs(a[i]); }
    }
    return minVal;
}

void arrMultScal(double res[], double a[], double b, int size){
    int i;
    for (i=0; i<size; i++) {
        res[i] = a[i] * b;
    }
    return;
}

void arrMult(double res[], double a[], double b[], int size){
    int i;
    for (i=0; i<size; i++) {
        res[i] = a[i] * b[i];
    }
    return;
}

double arrSum(double a[], int size){
    double res = a[0];
    int i;
    for (i=1; i<size; i++) {
        res += a[i];
    }
    return res;
}

void matMultCol(double res[], double *a[], double b[], int sizea, int sizeb){
    memset( res, 0, sizea * sizeof(double) );
    int j, k;
    for (j=0; j<sizea; j++){
        for (k=0; k<sizeb; k++){
            res[j] += a[k][j]*b[k];
        }
    }
    return;
}

void arrAdd(double res[], double a[], double b[], int size){
    int i;
    for (i=0; i<size; i++) {
        res[i] = a[i] + b[i];
    }
    return;
}

void arrSubtract(double res[], double a[], double b[], int size){
    int i;
    for (i=0; i<size; i++) {
        res[i] = a[i] - b[i];
    }
    return;
}

int anyNonNormal( double a[], int size){
    int res = 0;
    int i;
    for (i=0; i<size; i++) {
        if (isnan(a[i]) || isinf(a[i])){
            res = 1;
            return res;
        }
    }
    return res;
}

static int nbody_eq(
        double	dx[],
        double	x[],
        double   mus[],
        int      n
        ) {
    
    double temp[3] = { 0.0 };
    double rnorm3;
    int j;
    int k;
    int i;
    for(j = 0; j < n; j++){
        memset( temp, 0, sizeof(temp) );
        for(k = 0; k < n; k++){
            if (k==j){ continue; }
            rnorm3 = pow(NORM(x[k*3] - x[j*3], x[k*3+1] - x[j*3+1], x[k*3+2] - x[j*3+2]), 3);
            if (rnorm3 == 0.0){
                printf("Division by zero!\n");
		return -1;
            }
            for(i = 0; i < 3; i++){
                temp[i] += (x[k*3+i] - x[j*3+i])/rnorm3*mus[k];
            }
        }
        for(i = 0; i < 3; i++){
            dx[j*3+i] = temp[i];
        }
    }
    return 0;
}

static void cleanup_all(double **f, double *feval, double *fBphat, double *fBhat, double *fB, double *fBp, double *res1, double *res2){
    /* cleanup */
    int j;
    for (j = 0; j < 17; j++){
        free(f[j]);
    }
    free(feval);
    free(fBphat);
    free(fBhat);
    free(fB);
    free(fBp);
    free(res1);
    free(res2);
}

int nbodyRKN1210_c(double t_in[], double y_in[], double dy_in[], double mus[], int N, double y[], double dy[]){
    
    int j, k, sy0_m, res;
    sy0_m = N*3;
    int arrSize = sy0_m * sizeof (double);
        
    /* Initialize */
    double pw = 1.0/12.0;
    double abstol = 1e-14;
    double t0 = t_in[0];
    double tfinal = t_in[1];
    double t = t0;
    double hmax = fabs(tfinal - t0);
    double hmin = fabs(tfinal - t0)/1e12;
    memcpy( y, y_in,  arrSize );
    memcpy( dy, dy_in, arrSize  );
    int direction = (tfinal - t0)/fabs(tfinal - t0);
    double temp[17] = {0.0};
    
    double **f, *feval;
    feval = (double*) malloc(arrSize);
    memset( feval, 0, arrSize );
    
    f = (double**) malloc(17*sizeof(double*));
    for (j = 0; j < 17; j++){
        f[j] = (double*) malloc(arrSize);
        memset(f[j], 0, arrSize);
    }
    
    double *fBphat, *fBhat, *fB, *fBp;
    double *res1, *res2;
    double delta1, delta2, delta;
    
    fBphat = (double*) malloc(arrSize);
    memset( fBphat, 0, arrSize );
    
    fBhat = (double*) malloc(arrSize);
    memset( fBhat, 0, arrSize );
    
    fB = (double*) malloc(arrSize);
    memset( fB, 0, arrSize );
    
    fBp = (double*) malloc(arrSize);
    memset( fBp, 0, arrSize );
    
    res1 = (double*) malloc(arrSize);
    memset( res1, 0, arrSize );
    
    res2 = (double*) malloc(arrSize);
    memset( res2, 0, arrSize );
    
    /* initial step */
    res = nbody_eq(feval, y, mus, N);
    if (res != 0){
	cleanup_all(f,feval,fBphat,fBhat, fB, fBp, res1, res2);
	return -1;
    }
    memcpy(f[0], feval, arrSize );
    double h = pow(abstol, pw) /  MAX(MAX(arrAbsMax(dy, sy0_m), arrAbsMax(f[0], sy0_m)), 1e-4);
    h = MIN(hmax, MAX(h, hmin)) * direction;
    double h2 = h*h;

    double counter = 0.0;
    
    /* main loop */
    while (fabs(t-tfinal) > 0 ){
        
        /* correct if you're past final time */
        if (fabs(h) > fabs(tfinal-t)){
            h  = (tfinal - t);
            h2 = h*h;
        }
                    
        /* Compute the second-derivative */
        for (j = 0; j < 17; j++){
            matMultCol(res1, f, A[j], sy0_m, 17);
            arrMultScal(res1, res1, h2, sy0_m);
            arrMultScal(res2, dy, c[j]*h, sy0_m);
            arrAdd(res1, res1, res2, sy0_m);
            arrAdd(res1, res1, y, sy0_m);
            
            res = nbody_eq(feval, res1, mus, N);
	    if (res != 0){
		cleanup_all(f,feval,fBphat,fBhat, fB, fBp, res1, res2);
		return -1;
    	    }
            memcpy(f[j], feval, arrSize );
        }
    
        /* Estimate error */
        matMultCol(fBphat, f, Bphat, sy0_m, 17);
        matMultCol(fBhat, f, Bhat, sy0_m, 17);
        matMultCol(fBp, f, Bp, sy0_m, 17);
        matMultCol(fB, f, B, sy0_m, 17);
        
        arrSubtract(res1, fBhat, fB, sy0_m);
        arrMultScal(res1, res1, h2, sy0_m);
        delta1 = arrAbsMax(res1, sy0_m);
        
        arrSubtract(res2, fBphat, fBp, sy0_m);
        arrMultScal(res2, res2, h, sy0_m);
        delta2 = arrAbsMax(res2, sy0_m);
        
        delta = MAX(delta1, delta2)*fabs(h);
        
        /* update solution if error is acceptable */
        if (delta <= abstol){
            t = t + h;
            
            /*y  = y + h*dy + h^2*fBhat */  
            arrMultScal(res1, fBhat, h2, sy0_m);
            arrMultScal(res2, dy, h, sy0_m);
            arrAdd(res1, res1, res2, sy0_m);
            arrAdd(res1, res1, y, sy0_m);
            /* dy = dy + h*fBphat */ 
            arrMultScal(res2, fBphat, h, sy0_m);
            arrAdd(res2, res2, dy, sy0_m);
            
            /* Check for degeneracies */
            if ((anyNonNormal(y, sy0_m) != 0) || (anyNonNormal(dy, sy0_m) != 0)){
                printf("Failure at t = %6.6e:\n", t);
                printf("INF or NAN values encountered.");
                t = tfinal;
                cleanup_all(f,feval,fBphat,fBhat, fB, fBp, res1, res2);
		return -2;
            }else{
                memcpy(y, res1, arrSize );
                memcpy(dy, res2, arrSize );
            }
        }

        /* adjust step size */
        if (delta != 0){
            h = direction * MIN(hmax, 0.9*fabs(h)*pow(abstol/delta, pw));
            h2 = h * h;
            
            if (fabs(h) < hmin){
                printf("Failure at t = %6.6e:\nStep size fell below the minimum acceptible value of %6.6e \n", t, hmin);
                printf("Singularity Likely.\n");
                t = tfinal;
    		cleanup_all(f,feval,fBphat,fBhat, fB, fBp, res1, res2);
		return -3;
            }
        }
        counter++;
    }/* end main loop */
    
    cleanup_all(f,feval,fBphat,fBhat, fB, fBp, res1, res2);
    return 0; 
}

