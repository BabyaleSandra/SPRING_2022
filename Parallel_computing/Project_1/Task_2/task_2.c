/*                                             
Written by: Babyale Sandra                                                   
            Department of Mathematics                                           
            Boise State University                                                                      
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// 1) Generating random numbers
double random_number()
  {
    double f;
    f  = (double)rand() / (double)RAND_MAX ;
    return f;
  } 

void main(int argc, char* argv[]){

  long N = atol(argv[1]);

  //Declaring arrays
  double* x =  (double*) malloc(N*sizeof(double));
  double* x_bar =  (double*) malloc(N*sizeof(double));
  double* M  =  (double*) malloc(N*sizeof(double));

  //Declaring constants
  double sigma, x_xbar, x_bar_o = x[0]; 
  int i,k;
  
  //Intializing x_bar
  x_bar[0] = x_bar_o;

  //generating an array of random numbers
  for (i = 0; i <=N-1; ++i) x[i] = (double)random_number();
  

  // Welford's algorithm                       
   for (k = 1; k<=N-1; ++k)
   {
      x_xbar = (x[k] - x_bar[k-1]);
      x_bar[k] = x_bar[k-1] + x_xbar/k;
      M[k] = M[k-1] + (x[k] - x_bar[k])*x_xbar;
   }

  //compute the standard deviation                                              
  sigma = sqrt(M[N-1]/N);

  printf("std_dev = %e\n",sigma);

}