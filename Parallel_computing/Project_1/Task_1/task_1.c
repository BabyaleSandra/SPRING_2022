#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// 1) Generating random numbers
double random_number()
  {
    double f;
    f  = (double)rand() / (double)RAND_MAX ;
    return f;
  } 


 
 
void main(int argc, char* argv[]){
  
long N = atol(argv[1]);
 MPI_Init(&argc, &argv);
 int i, nproc, irank, root_rank=0;
 

 double* x_loc;
 double* x=  (double*) malloc(N*sizeof(double));
 
  
     
 //Getting array of random numbers

 
  
 //find out communicator size and rank number
 MPI_Comm_size(MPI_COMM_WORLD,&nproc);
 MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  
 //start timer
 double start_timer = MPI_Wtime();
 

 // 2) LET ALL RANKS KNOW HOW MANY POINTS THERE ARE 
  
  double start_timer_broadcast = MPI_Wtime();
   

  MPI_Bcast(&N, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
  double broadcast_comm_time = MPI_Wtime() -  start_timer_broadcast ;
  
  double start_comput_timer1  = MPI_Wtime();
  
  // 3) Compute the number of points per rank as N_loc and generating x_loc
  int N_loc = N/nproc;
  x_loc = (double*) malloc(N_loc*sizeof(double));
  
 // 4) Creating array of random numbers for each rank
 for(i=0; i<N_loc; i++){
   x_loc[i] = (double)random_number();
    
 }
 
 double end_comput_timer1  = MPI_Wtime() - start_comput_timer1 ;
  
  
 
  
 // 5) --------- COMPUTE THE SUM OF DATA POINTS --------------------------
  
  double start_comput_timer2 = MPI_Wtime();
  double mu_loc = 0;
  for(i=0;i<N_loc;i++)
  {
    mu_loc += x_loc[i];
  }
 double end_comput_timer2 =  MPI_Wtime()- start_comput_timer2;
  
  // 6) --------- COMPUTE THE GLOBAL SUM ----------------------------------

  double start_timer_allreduce = MPI_Wtime();
  double mu;
  MPI_Allreduce(&mu_loc, &mu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  double allreduce_comm_time = MPI_Wtime() - start_timer_allreduce;
  // 7) --------- COMPUTE THE MEAN VALUE ----------------------------------
 double start_comput_timer3  = MPI_Wtime();
  mu = mu/N;
  // 8) --------- COMPUTE THE SUM OF DEVIATIONS ---------------------------
  //remember to use an appropriate loop size and variable name
  double sigma_loc=0;
  for(i=0;i<N_loc;i++){
   sigma_loc += (x_loc[i] - mu)*(x_loc[i] - mu);
  }
  double end_comput_timer3 =  MPI_Wtime() - start_comput_timer3;

  
  // 9) --------- COMPUTE GLOBAL SUM OF SIGMA_LOC -------------------------
  double start_timer_reduce = MPI_Wtime();
  double sigma;
  MPI_Reduce(&sigma_loc, &sigma, 1, MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  double reduce_comm_time = MPI_Wtime() - start_timer_reduce;
  // 10) --------- COMPUTE THE STANDARD DEVIATIONS -------------------------
  //only on rank 0
  sigma = sqrt(sigma/N);

  //stop timer for computation
  double stop_timer = MPI_Wtime();
  double elapsed_time = stop_timer - start_timer;
  
  // Getting the communication time and computational time
  double communication_time = broadcast_comm_time + allreduce_comm_time + reduce_comm_time;
  double computational_time = end_comput_timer3+ end_comput_timer2+ end_comput_timer1; 
  
  double timers[9] = {broadcast_comm_time, allreduce_comm_time,reduce_comm_time,communication_time,end_comput_timer1,end_comput_timer2,end_comput_timer3,computational_time,elapsed_time};

  double elapsed_time_global[9];

  MPI_Reduce(timers, elapsed_time_global, 9, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

  if(irank==root_rank){
    printf("%d,%ld,%f,%f",nproc,sigma);
    for(i=0; i<9;i++){
      printf(",%f",timers[i] );
    }
    printf("\n");
  }

  MPI_Finalize();
}
