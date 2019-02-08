#include <math.h>
#include <string.h>
#include "nbodyC.h"

/*N-body equations.
Calculate the mutual accelerations of a collection of N bodies defined by 
gravitational parameters (G*M) stored in self.mu

Args:
    dx (double array):
        Return values (x'')
    x (double array)
        State value (body positions 3nx1)
    mus (double array)
        Gravitational parameters (nx1)
    n (scalar)
        Number of bodies
*/

int nbodyC(double dx[],double x[],double mus[],int n)
{    
    double temp[3] = { 0.0 };
    double rnorm3;
    int j;
    int k;
    int i;
    for(j = 0; j < n; j++){
        memset( temp, 0, sizeof(temp) );
        for(k = 0; k < n; k++){
            if (k==j){ continue; }
	    rnorm3 = pow(pow(x[k*3] - x[j*3],2.) + pow(x[k*3+1] - x[j*3+1],2.) + pow(x[k*3+2] - x[j*3+2],2.),3./2.);
            if (rnorm3 == 0.0){
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
