#include <R.h>



void fullconvolvePaired(double *pb, int *n, double *result){
  
  int L = 2; // length of previous interation
  double p,q;
  double tmp_head,tmp_tail;
  
  result[0] = 1-pb[0];
  result[1] = pb[0];
  
  for(int i = 1; i < *n; i++){
    p = pb[i];
    q = 1 - pb[i];
    
    //Edge cases
    result[L] = p * result[L-1]; 
    tmp_head = result[0];
    result[0] = tmp_head * q;
    
    //Remaining cases
    for(int j = 1; j < L; j++){
      tmp_tail = result[j];
      result[j] = p * tmp_head + q * result[j];
      tmp_head = tmp_tail;
    }
    
    L++;
  }
  
  return;
  
}

