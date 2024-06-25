/*
The code solves a 2D Poisson equation in a square domain using finite difference method

Written by: Michal A. Kopera
            Department of Mathematics
	    Boise State University
	    03/06/2021
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

// *** Function interfaces ***

// getArgs_mpi reads input parameters Nfrom command line
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm);
  
//compute integral L2 error
double L2_error(int N_loc, int N, double dx, double L, double* data, double* reference);

// write results to file
void write_results(double *u, int N, int N_loc, int nproc, int irank, double h);


// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, id_ghost;
  int id_lft, id_rgt, id_top, id_bot;

  int irank, nproc;
  
  double eps = 1e-8; //error tolerance for Jacobi method
  
  //declare variables for tracking requests
  MPI_Request request_send_rgt;
  MPI_Request request_send_lft;
  MPI_Request request_send_bot;
  MPI_Request request_send_top;

  MPI_Request request_recv_rgt;
  MPI_Request request_recv_lft;
  MPI_Request request_recv_bot;
  MPI_Request request_recv_top; 
  
  //initialize MPI
  MPI_Init(&argc,&argv); 
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  //compute row and column index in the rank grid
  int q = (int)sqrt(nproc); //assume nproc is a square number
  int rank_row = irank/q;
  int rank_col = irank%q;
  
  //read command line arguments
  getArgs(&N, argc, argv, irank, MPI_COMM_WORLD);

  // define domain size (assume both x and y directions are the same)
  double L = 1.0;
  
  //compute the grid spacing
  double h = L/(N-1);

  //compute local problem size
  double N_loc = N/q;
  
  //allocate arrays
  double *u       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *u_new   = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *u_exact = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *f       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *x       = malloc((N_loc+2)*sizeof(double));
  double *y       = malloc((N_loc+2)*sizeof(double));

  //allocate buffers
  double *send_buffer_rgt = malloc(N_loc*sizeof(double));
  double *send_buffer_lft = malloc(N_loc*sizeof(double));
  double *send_buffer_bot = malloc(N_loc*sizeof(double));
  double *send_buffer_top = malloc(N_loc*sizeof(double));

  double *recv_buffer_rgt = malloc(N_loc*sizeof(double));
  double *recv_buffer_lft = malloc(N_loc*sizeof(double));
  double *recv_buffer_bot = malloc(N_loc*sizeof(double));
  double *recv_buffer_top = malloc(N_loc*sizeof(double));

  
  //initialize x and y
  for(i=1;i<N_loc+1;i++){
    x[i] = (i-1)*h + rank_col*N_loc*h;
    y[i] = (i-1)*h + rank_row*N_loc*h;
  }
  
  //initialize array u
  for(j=1;j<N_loc+1;j++){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(i,j,N_loc);
	u_exact[id] = sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);

	f[id] = -8*M_PI*M_PI*sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);

	u[id] = 0.0;
      }
  }

  //initialize BC in ghost cells

  if(rank_row==0){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,0,N_loc); //bottom BC
      u[id] = 0.0;
    }
  }

  if(rank_row==q-1){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,N_loc+1,N_loc); //top BC
      u[id] = 0.0;
    }
  }

  if(rank_col==0){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(0,i,N_loc); //left BC
      u[id] = 0.0;
    }
  }

  if(rank_col==q-1){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(N_loc+1,i,N_loc); //right BC
      u[id] = 0.0;
    }
  }
       
  //compute initial error
  double error = L2_error(N_loc, N, h, L, u, u_exact);

  int iter=0;

  double time_start = MPI_Wtime();

  double time_comm_start, time_comp_start, time_edge_start, time_l2_start;

  double time_comm = 0;
  double time_comp = 0;
  double time_edge = 0;
  double time_l2   = 0;
  
  //begin Jacobi iteration
  while(error>eps){

    //initialize communication

    time_comm_start = MPI_Wtime();
    //pack data to send buffers
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(0,i,N_loc);
      send_buffer_lft[i] = u[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(N_loc+1,i,N_loc);
      send_buffer_rgt[i] = u[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,0,N_loc);
      send_buffer_bot[i] = u[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,N_loc+1,N_loc);
      send_buffer_top[i] = u[id];
    }

    //compute neighbor indices
    int rank_lft  = irank-1;
    int rank_rgt  = irank+1;
    int rank_bot  = irank-q;
    int rank_top  = irank+q;
    
    //initialize receives
    if(rank_col>0)   MPI_Irecv(recv_buffer_lft, N_loc, MPI_DOUBLE, rank_lft, 101, MPI_COMM_WORLD, &request_recv_lft);
    if(rank_col<q-1) MPI_Irecv(recv_buffer_rgt, N_loc, MPI_DOUBLE, rank_rgt, 102, MPI_COMM_WORLD, &request_recv_rgt);
    if(rank_row>0)   MPI_Irecv(recv_buffer_bot, N_loc, MPI_DOUBLE, rank_bot, 103, MPI_COMM_WORLD, &request_recv_bot);
    if(rank_row<q-1) MPI_Irecv(recv_buffer_top, N_loc, MPI_DOUBLE, rank_top, 104, MPI_COMM_WORLD, &request_recv_top);

    //initialize sends
    if(rank_col<q-1) MPI_Isend(send_buffer_rgt, N_loc, MPI_DOUBLE, rank_rgt, 101, MPI_COMM_WORLD, &request_send_rgt);
    if(rank_col>0)   MPI_Isend(send_buffer_lft, N_loc, MPI_DOUBLE, rank_lft, 102, MPI_COMM_WORLD, &request_send_lft);
    if(rank_row<q-1) MPI_Isend(send_buffer_top, N_loc, MPI_DOUBLE, rank_top, 103, MPI_COMM_WORLD, &request_send_top);
    if(rank_row>0)   MPI_Isend(send_buffer_bot, N_loc, MPI_DOUBLE, rank_bot, 104, MPI_COMM_WORLD, &request_send_bot);


    time_comp_start = MPI_Wtime();
    //go over interior points 2:N_loc-1
    for(j=2;j<N_loc-1;j++){
      for(i=2;i<N_loc-1;i++){
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top

	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }


    
    //enforce Dirichlet BC by copying values from ghosts
    if(rank_row ==0){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(i,1,N_loc); //bottom boundary
	id_ghost = ID_2D(i,0,N_loc); //bottom ghost
	u_new[id] = u[id_ghost];
      }
    }

    if(rank_row == q-1){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(i,N_loc,N_loc); //top boundary
	id_ghost = ID_2D(i,N_loc+1,N_loc); //top ghost
	u_new[id] = u[id_ghost];
      }
    }

    if(rank_col == 0){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(1,i,N_loc); //left boundary
	id_ghost = ID_2D(0,i,N_loc); //left ghost
	u_new[id] = u[id_ghost];    
      }
    }

    if(rank_col==q-1){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(N_loc,i,N_loc); //right boundary
	id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost
	u_new[id] = u[id_ghost];
      }
    }
    time_comp = time_comp + (MPI_Wtime() - time_comp_start);
    
    //check if communication complete
    if(rank_col>0){
      MPI_Wait(&request_recv_lft,MPI_STATUS_IGNORE); //wait until left recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(0,i,N_loc); //left ghost
	u[id_ghost] = recv_buffer_lft[i]; //load receive buffer to left ghost
      }      
    }

    if(rank_col<q-1){
      MPI_Wait(&request_recv_rgt,MPI_STATUS_IGNORE); //wait until right recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost
	u[id_ghost] = recv_buffer_rgt[i]; //load receive buffer to right ghost
      }      
    }

    if(rank_row>0){
      MPI_Wait(&request_recv_bot,MPI_STATUS_IGNORE); //wait until bottom recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,0,N_loc); //bottom ghost
	u[id_ghost] = recv_buffer_bot[i]; //load receive buffer to bottom ghost
      }      
    }
    
    if(rank_row<q-1){
      MPI_Wait(&request_recv_top,MPI_STATUS_IGNORE); //wait until top recv is complete
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,N_loc+1,N_loc); //top ghost
	u[id_ghost] = recv_buffer_top[i]; //load receive buffer to top ghost
      }      
    }

    time_comm = time_comm + (MPI_Wtime() - time_comm_start);
    //complete computation at the rank edges (except the domain boundary conditions)

    time_edge_start = MPI_Wtime();
    if(rank_col>0){ //compute left edge
      for(j=1;j<N_loc+1;j++){
	i=1;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }

    if(rank_col<q-1){ //compute right edge
      for(j=1;j<N_loc+1;j++){
	i=N_loc;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }


    if(rank_row>0){ //compute bottom edge
      for(i=1;i<N_loc+1;i++){
	j=1;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }

    
    if(rank_row<q-1){ //compute top edge
      for(i=1;i<N_loc+1;i++){
	j=N_loc;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }
    time_edge = time_edge + (MPI_Wtime() - time_edge_start);
    // computation done here

    
    //compute error
    time_l2_start = MPI_Wtime();
    error = L2_error(N_loc,N, h, L, u_new, u);
    time_l2 = time_l2 + (MPI_Wtime() - time_l2_start);
    iter++;
    
    //update solution 
    for(j=1;j<N_loc+1;j++){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(i,j,N_loc);
	u[id] = u_new[id];
      }
    }
        
  }

    double time_elapsed = MPI_Wtime() - time_start;
  
    if(irank==0) printf("N = %d, iter = %d, eps = %e, error = %e, time_elapsed = %e, time_comp = %e, time_comm = %e, time_edge = %e, time_l2 = %e\n",N,iter,eps,error,time_elapsed, time_comp, time_comm, time_edge, time_l2);

    //write_results(u, N, N_loc,nproc, irank , h);
  
    
#ifdef DEBUG // print the content of array u including ghosts
  for(j=0; j<N+2; j++){
    for(i=0; i<N+2; i++){
      id = ID_2D(i,j,N);
      printf("%f\t",u[id]);
    }
    printf("\n");
  }
#endif
  
  free(u);
  free(u_new);
  free(u_exact);
  free(f);
  free(x);
  free(y);

  MPI_Finalize();
  return 0;
}


// *** Function definitions ***
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm)
{

  if(irank==0){
    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./poisson N \n");
	MPI_Finalize ( );
	exit ( 1 ); // exit with error code 1
      }
    else
      {
	//Use input argument
	*N = atoi(argv[1]);
      }
  }
  MPI_Bcast(N,1,MPI_INT,0,comm);
}


double L2_error(int N_loc, int N, double dx, double L, double* data, double* reference){

  double error=0.0;
  double error_total = 0.0;
  int id;
  
  for(int j=1;j<N_loc+1;j++){
    for(int i=1;i<N_loc+1;i++){
      id = ID_2D(i,j,N_loc);
      error += pow((data[id]-reference[id]),2);
    }
  }
  MPI_Allreduce(&error,&error_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return sqrt(error_total)/N;
}

void write_results(double *u, int N, int N_loc, int nproc, int irank, double h){

  int i,j, id;

  double *u_local       = malloc((N_loc)*(N_loc)*sizeof(double));
  double *u_global      = malloc((N)*(N)*sizeof(double));
  double *u_write      = malloc((N)*(N)*sizeof(double));

  double x,y;
  
  //pack the data for gather (to avoid sending ghosts)
  int id_loc = 0;
  for(j=1;j<N_loc+1;j++){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,j,N_loc);
      u_local[id_loc] = u[id];
      id_loc++;
    }
  }
  
  MPI_Gather(u_local,id_loc,MPI_DOUBLE,u_global,id_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);

  //unpack data so that it is in nice array format

  int id_write, id_global;
  int p_row, p_col;
  int q = sqrt(nproc);

  if(irank==0){
    for(int p=0; p<nproc;p++){
      p_row = p/q;
      p_col = p%q;
      for(j=0;j<N_loc;j++){
	for(i=0;i<N_loc;i++){
	  id_global = p*N_loc*N_loc + j*N_loc + i;
	  id_write  = p_row*N_loc*N_loc*q + j*N_loc*q + p_col*N_loc + i;
	  u_write[id_write] = u_global[id_global];
	}
      }
    }

#ifdef DEBUG
    //writing to screen
    printf("Checking unpacked array again\n");
    for(j=0;j<N;j++){
      for(i=0;i<N;i++){
	id = j*N+i;
	printf("%f ",u_write[id]);
      }
      printf("\n");
    }
#endif

    
    //write to file
    FILE *f;
    f = fopen("results_poisson.dat","w"); //open file
    fprintf(f,"x, y, u\n");

    for(j=0; j<N; j++){
      for(i=0; i<N; i++){
	id = j*N + i;
	x = i*h; //I am a bit lazy here with not gathering x and y
	y = j*h;
	fprintf(f,"%f, %f, %f\n",x,y,u_write[id]);
      }
    }
    fclose(f);
  }
  free(u_local); free(u_global); free(u_write);
  
}
