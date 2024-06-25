/*                                             
Written by: Babyale Sandra                                                   
            Department of Mathematics                                           
            Boise State University                                                                      
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>  

// 1) GENERATING RANDOM NUMBERS
double random_number()
{
    double f;
    f  = (double)rand() / (double)RAND_MAX ;
    return f;
} 

void main(int argc, char* argv[]){
  
  MPI_Init(&argc, &argv);

  // Declaring constants
  int  N_loc, nproc, irank, root_rank = 0, i, k;
  double sigma, x_xbar,M_global, sigma1, x_bar_g, x_bar_o, x_bar_k,xbar_o;
  long N = atol(argv[1]);

  
  // find out communicator size and rank number                                              
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);                                                     
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  //start timer
 double start_timer = MPI_Wtime();
 
 double start_comput_timer1  = MPI_Wtime();
 // 2) COMPUTE THE NUMBER OF POINTS PER RANK - N_LOC 
   N_loc = N/nproc;

   // Declaring arrays
  double* x =  (double*) malloc(N_loc*sizeof(double));
  double* x_bar =  (double*) malloc(N_loc*sizeof(double));
  double* M  =  (double*) malloc(N_loc*sizeof(double));
  double* x_bar_global =(double*) malloc(nproc*sizeof(double));

 
 // 3) GENERATING AN ARRAY OF RANDOM NUMBERS
  srand(irank);
  for (i = 0; i < N_loc; i++) x[i] = (double)random_number();
  
 // Initializing x_bar   
  x_bar_o = x[0];
  x_bar[0] = x_bar_o; 

 // 4) WELFORD'S ALGORITHM 
for (k = 1; k<N_loc; k++)
  {
      x_xbar = (x[k] - x_bar[k-1]);
      x_bar[k] = x_bar[k-1] + x_xbar/k;
      M[k] = M[k-1] + (x[k] - x_bar[k])*x_xbar;
  }
  x_bar_k = x_bar[N_loc-1];
  double end_comput_timer1  = MPI_Wtime() - start_comput_timer1 ;

 double start_timer_gather = MPI_Wtime();
  MPI_Gather(&x_bar_k,1,MPI_DOUBLE,x_bar_global,1,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);
 double gather_comm_time = MPI_Wtime() -  start_timer_gather ;
  
  double start_timer_reduce = MPI_Wtime();
  MPI_Reduce(&M[N_loc-1],&M_global,1,MPI_DOUBLE,MPI_SUM,root_rank,MPI_COMM_WORLD);
  double reduce_comm_time = MPI_Wtime() -  start_timer_reduce ;

  double start_comput_timer2 = MPI_Wtime();
  x_bar_o = x_bar_global[0];
  for (i = 1; i<nproc; i++)
    {
       x_bar_g = x_bar_o + (- x_bar_o + x_bar_global[i])*(i*N_loc)/((i+1)*N_loc);
       M_global += (pow((x_bar_o - x_bar_global[i]),2))*(i*pow(N_loc,2)/((i+1)*N_loc));
       x_bar_o = x_bar_g;
    }
  
  //compute the standard deviation 
  sigma = sqrt(M_global/N);
  double end_comput_timer2  = MPI_Wtime() - start_comput_timer2; 

  //stop timer for computation
  double stop_timer = MPI_Wtime();
  double elapsed_time = stop_timer - start_timer;

  // Getting the communication time and computational time
  double communication_time = gather_comm_time +  reduce_comm_time;
  double computational_time = end_comput_timer2+ end_comput_timer1; 
  
  double timers[7] = {gather_comm_time,reduce_comm_time,communication_time,end_comput_timer1,end_comput_timer2,computational_time,elapsed_time};
  double elapsed_time_global[7];

  MPI_Reduce(timers, elapsed_time_global, 5, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

  if(irank==root_rank){
    printf("%d,%f",nproc,sigma);
    for(i=0; i<7;i++){
      printf(",%f",timers[i] );
    }
    printf("\n");
  }

  MPI_Finalize();
}
   
  //double timers[7] = {gather_comm_time,reduce_comm_time,communication_time,end_comput_timer1,end_comput_timer2,computational_time,elapsed_time};
  //double elapsed_time_global[7];