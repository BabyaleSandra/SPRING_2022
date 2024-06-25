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

  // compute how many points each rank has
  long N_loc = N/nproc;

  //allocate arrays
  double* u    = (double*) malloc(N_loc*sizeof(double));
  
  // compute the interval size dx - we integrate from 0 to 1
  double dx = 1.0/(N-1);

  // store the value of exact solution
  double pi_exact = M_PI;

  //start timer for initialization
  double start_timer_init = MPI_Wtime();
  
  //create initial condition
  init_u(u,dx,N_loc,irank);

  //stop timer for initialization
  double stop_timer_init = MPI_Wtime();
    
  //start timer for computation
  double start_timer_comput = MPI_Wtime();

  //call finite difference function
  double I_loc = trapez(u,N_loc,dx);

  //stop timer for computation
  double stop_timer_comput = MPI_Wtime();

  //compute the elapsed time
  double elapsed_time_init = stop_timer_init - start_timer_init;
  double elapsed_time_comput = stop_timer_comput - start_timer_comput;
  
  //print the elapsed time
  //printf("Problem size N = %ld, N_loc = %ld, irank = %d, I_loc = %e,  Elapsed time: initialization = %e s, computation = %e s\n",N,N_loc,irank,I_loc,elapsed_time_init, elapsed_time_comput);

  
  //which rank will be the root
  int root_rank=0;
  
  double I; 
  MPI_Reduce(&I_loc, &I, 1, MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  
  double elapsed_time_comput_global, elapsed_time_init_global;
  MPI_Reduce(&elapsed_time_comput, &elapsed_time_comput_global, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(&elapsed_time_init,   &elapsed_time_init_global,   1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

  if(irank==root_rank){
    //printf("Problem size N = %ld, irank = %d, I = %16.14f,  Elapsed time: initialization = %e s, computation = %e s\n",N,irank,I,elapsed_time_init_global, elapsed_time_comput_global);
  printf("%d,%e\n",nproc,elapsed_time_comput_global);
  
  };



  

  //free allocated memory
  free(u);

  MPI_Finalize();
}
