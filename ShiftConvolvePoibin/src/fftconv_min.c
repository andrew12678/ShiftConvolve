#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <fftw3.h>
#include "minfft.h"
#include <stdlib.h>
#include <string.h>
#include <R.h>

/*
 * This is the C implemented 'streamlined' approach of convoluting a given probability vector in the Fourier domain
*/
void fftconvPairs(double *p, int *n, double complex *A, double complex *B, double complex *U,  double *result, double *in, double complex *w_vector, double complex *out){
    int H = 4; //Each column belongs to 1 bernoulli prob and starts off with 4 values (p, 1-p, 0, 0)
    int L = *n; //Number of columns (different bernoulli probabilities to start with)

    minfft_aux *a; //Necessary auxillary data for minFFT package
    a = minfft_mkaux_realdft_1d(H);

    in[2] = 0; //We are always padding the last two entries of the columns with zeroes
    in[3] = 0;
    for(int i = 0; i < L; i++){
      in[0] = p[i];
      in[1] = p[i+L];
      minfft_realdft(in, out, a);
      A[i*H] = out[0]; //Storing the Fourier domain transformations (complex numbers) in complex array A
      A[i*H+1] = out[1];
      A[i*H+2] = out[2];
      A[i*H+3] = conj(out[1]);
    }

    minfft_free_aux(a);

    double complex w = 0+1*I; //The imaginary number i
    minfft_aux *a2, *a3;
    a2 = minfft_mkaux_realdft_1d(H);
    a3 = minfft_mkaux_dft_1d(H);
    while(L!=2){ //Iterate until there are only two columns remaining
      int c; //New number of columns (ceiling)
      int d; //New number of columns (floor)
      if(L % 2 == 0){ //If L even, c and d should be the same
        c = L/2;
        d = L/2;
      }else{ //If L odd, c is the ceiling, d is the floor
        c = L/2 + 1;
        d = L/2;
      }

      for(int i = 0; i < d; i++){
        for(int j = 0; j < H; j++){
          //B is a complex array declared with sufficient memory to handle all potential new dimensions (c,d)
          //and is designated to store the pointwise multiplication of paired columns
          //The memory for B is reused in each iteration.
          B[i*H+j] = A[2*i*H + j] * A[(2*i+1)*H +j];
        }
      }

      if(d != c){ //If L odd, logic set up by c and d earlier
        for(int i = 0; i < H; i++){
          B[(c-1)*H+i] = A[(L-1)*H+i];
        }
      }

      if(H != 4){ //Dynamically resize memory of our fftw arrays when the number of entries in each column changes at end of loop
        minfft_free_aux(a2);
        a2 = minfft_mkaux_realdft_1d(H); //Auxiliary Data required for Inverse Complex -> Real FFT
        minfft_free_aux(a3);
        a3 = minfft_mkaux_dft_1d(H); //Auxiliary Data required for Complex -> Complex FFT

      }

      //Pre-compute omega
      for(int j = 0; j < H; j++){
        w_vector[j] = cexp(-w*M_PI/H*j);
      }

      // Inverse FFT to retreive pmf values into U array
      for(int i = 0; i < c; i++){
        minfft_invrealdft(B + i*H, in, a2); //Inverse Complex -> Real FFT
        //Executing pointwise multiplication in the same loop for efficiency
        for(int j = 0; j < H; j++){
          U[i*H + j] = in[j] / H * w_vector[j];
        }


      }

      //Pointwise multiplication of U*v and FFT
      for(int i = 0; i < c; i++){
        //Complex --> Complex 1D-FFT of the pointwise multiplication
        minfft_dft(U + i*H, out, a3); //Complex -> Complex FFT
        memcpy(U + i*H, out, sizeof(minfft_cmpl)*H);

      }

      //Advancing odd and even indicies
      for(int i = 0; i < 2*c*H; i++){
        if(i % 2 == 0){
          A[i] = B[i/2];
        }else{
          A[i] = U[i/2];
        }
      }

      //Updating size of columns and number of columns
      H *= 2;
      L = c;
    }


    //Final multiplication and inverse FFT to retreive pmf

    for(int i = 0; i < H; i++){
      A[i] = A[i]*A[H+i];
    }

    minfft_free_aux(a2);
    minfft_free_aux(a3);
    minfft_aux *a4;
    a4=minfft_mkaux_dft_1d(H);
    minfft_invdft(A,out,a4); //Inverse complex -> complex FFT
    for(int i = 0; i < H; i++){
      result[i] = cabs((out[i]/H));
    }

    minfft_free_aux(a4);
    return;
}


/*
int main(){
    double p[] = {0.05631382, 0.61668317, 0.61238685, 0.88876782, 0.63543906, 0.94368618, 0.38331683, 0.38761315, 0.11123218, 0.36456094};
    int n = 5;
    double complex A[20];
    double complex B[20];
    double complex U[20];
    double result[20];
    fftconvPairs(p, &n, A, B, U, result);

}
*/



