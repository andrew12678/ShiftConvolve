#include <math.h>
#include <stdlib.h>
#include <R.h>

double logsum(double t0, double p1, double p2, double *xmax){
    //double xmax;
    if(p1 >= p2+t0){
        *xmax = p1;
    }else{
        *xmax = p2 + t0;
    }

    return *xmax + log(exp(p1 - *xmax) + exp(p2 + t0 - *xmax));
}

void computeMGF(double *p, int *n, double *t0, double *mgf, double *shift, double *xmax){

    for(int i = 0; i < *n; i++){
        mgf[i] = logsum(*t0,p[i], p[i+*n], xmax);
        shift[i] = p[i] - mgf[i];
        shift[i+*n] = p[i+*n] + (*t0) - mgf[i];

    }

    return;


}
