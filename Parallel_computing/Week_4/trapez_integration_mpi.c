/***************************
trapez_integration computes an approximation of a definite integral of a function
given by a set of discrete points using trapezoidal rule

I = sum(0.5*dx*(u(i+1)+u(i)))

Written by: Michal A. Kopera
            Department of Mathematics
            Boise State University
            1/25/2021

 **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

/*
trapez comoutes the approximation to the definite integral I to function u given by a set of discrete points

Inputs:
u  - an array of values of function u at discrete points
N  - number of points in the u array
dx - distance between points (assunes equidistant distribution)

Outputs:
I  - the value of the integral (returned through function name)
*/
double trapez(double* u, long N, double dx)
{
  long i; // local index for traversing arrays

  //initialize the integral variable I
  double I = 0;

  // go over interior points and compute second-order finite difference
  for (i=0; i<N-1; ++i)
    {
      I =  I + 0.5*dx*(u[i+1] + u[i]);
    }

  return I;
}

/*
init_u initializes array u with some data
Here we use a 4.0/(1+x*x) function

Inputs:
N  - the number of points
dx - the distance between points
 
Outputs:
u  - the values of u at specified points
*/

void init_u(double* u, double dx, long N_loc, int irank)
{
  long i, i_loc;
  double x;
  
  //for each point compute the value of our function at i*dx
  for (i_loc=0; i_loc<N_loc; ++i_loc)
    {
      i = N_loc*irank + i_loc;
      x = i*dx;
      u[i_loc] = 4.0/(1+x*x);
    }
}

/* main function 
Creates the function which derivative is computed and measures the execution time of the finite difference computation. 

Assumes that it receives N as the first (after the executable name) command line argument
Inputs:
N - number og points used to approximate function u

*/

int main(int argc, char* argv[])
{

  
  //Initialize MPI - create all available processes
  MPI_Init(&argc, &argv);

  //Variables to store communicator size and rank number
  int nproc, irank; 

  //Check communicator size
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  // get the number of points from input parameters
  long N = atol(argv[1]);

  // **** compute how many points each rank has
  long N_loc = (N-1)/nproc +1;
  printf("N_loc = %d, irank = %d \n",N_loc,irank);
 
  //allocate arrays
  double* u    = (double*) malloc(N_loc*sizeof(double));
  
  // compute the interval size dx - we integrate from 0 to 1
  double dx = 1.0/(N-1);

  // store the value of exact solution
  double pi_exact = M_PI;
  
  //create initial data
  init_u(u,dx,N_loc,irank);

  // **** compute integral
  double I = trapez(u,N_loc,dx);

  // **** Add partial resuults and store on rank 0
    

  //compute the error
  double pi_error = fabs(pi_exact-I);

  
    printf("Problem size N = %ld, I = %e, error = %e\n",N,I, pi_error);
  

  

  //free allocated memory
  free(u);

  MPI_Finalize();
}
