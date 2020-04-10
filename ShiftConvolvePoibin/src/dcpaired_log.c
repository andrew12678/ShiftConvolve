#include <R.h>
#include <math.h>

/*
double logsum2(double p1, double p2){
  double xmax;
  if(p1 >= p2){
    xmax = p1;
  }else{
    xmax = p2;
  }
  return xmax + log(exp(p1 - xmax) + exp(p2 - xmax));
}
*/

void fullconvolvePairedLog(double *pb, int *n, double *result){
  
  int L = 2; // length of previous interation
  double p,q;
  double tmp_head,tmp_tail;
  double a,b;
  
  result[0] = log(1-pb[0]);
  result[1] = log(pb[0]);
  
  for(int i = 1; i < *n; i++){
    p = log(pb[i]);
    q = log(1 - pb[i]);
    
    //Edge cases
    result[L] = p + result[L-1]; 
    tmp_head = result[0];
    result[0] = tmp_head + q;
    
    //Remaining cases
    for(int j = 1; j < L; j++){
      tmp_tail = result[j];
      b = q + tmp_tail;
      a = p + tmp_head;
      if(a>b){
        result[j] = a + log(1 + exp(b - a));
      }else{
        result[j] = b + log(1 + exp(a - b));
      }
      //result[j] = logsum2(p + tmp_head, q + result[j]);
      tmp_head = tmp_tail;
    }
    
    L++;
  }
  
  return;
  
}
