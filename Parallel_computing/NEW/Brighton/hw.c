#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NX  18 // This is the number of x cells in each process
#define NY  18 // This is the number of y cells in each processes
#define  R  16 //# of process
int nsx = 4;
int nsy = 4;
int block_ID[2];

// local Variables
double qu[NY][NX],qu_old[NY][NX],qd[NY][NX],qd_old[NY][NX],
     qr[NY][NX],qr_old[NY][NX],ql[NY][NX],ql_old[NY][NX];

double lu[NY][NX],lu_old[NY][NX],ld[NY][NX],ld_old[NY][NX],
       lr[NY][NX],lr_old[NY][NX],ll[NY][NX],ll_old[NY][NX];

double bd_condl[NY][2], bd_condr[NY][2], bd_condd[NX][2], bd_condu[NX][2];

double p[NY][NX],p_old[NY][NX];

double betau[NY][NX],betad[NY][NX],betar[NY][NX],betal[NY][NX];
double Ku[NY][NX],Kd[NY][NX],Kl[NY][NX],Kr[NY][NX];
double f[NY][NX];

double perm[NY][NX];

double h, maxiter, cfactor, size, tol; // note the grids are squares

// List all functions here

void set_subdomain_ID(int block_ID[2], int rank, int nsx) {
  block_ID[0] = rank%nsx;
  block_ID[1] = rank/nsx;
}

void set_domain_type(int block_ID[2], int block_type[2], int nsx, int nsy) {
 /* 4 corners */
  if(block_ID[0] == 0 && block_ID[1] == 0) {
    block_type[1] = 0;
  } else if(block_ID[0] == nsx-1 && block_ID[1] == 0) {
    block_type[1] = 2;
  } else if(block_ID[0] == nsx-1 && block_ID[1] == nsy-1) {
    block_type[1] = 8;
  } else if(block_ID[0] == 0 && block_ID[1] == nsy-1) {
      block_type[1] = 6;
  } else if (block_ID[1] == 0) {
    block_type[1] = 1;
  } else if(block_ID[0] == nsx-1) {
    block_type[1] = 5;
  } else if(block_ID[1] == nsy-1) {
      block_type[1] = 7;
  } else if(block_ID[0] == 0) {
    block_type[1] = 3;
  } else {
    block_type[1] = 4;
  }
}

void set_betas(double Kl[NY][NX], double Kr[NY][NX], double Kd[NY][NX], double Ku[NY][NX], double perm[NY][NX],
  double betal[NY][NX], double betar[NY][NX], double betad[NY][NX], double betau[NY][NX], double cfactor, double h) {
  for(int i=0;i<NY;i++) //initializing K-alpha field
  {
    for(int j=0;j<NX;j++)
    {
      if(i==0 || i==NY || j==0 || j==NX)
      {
        Ku[i][j]=0;
        Kd[i][j]=0;
        Kl[i][j]=0;
        Kr[i][j]=0;
      }
      else
      {
        Ku[i][j]=2*perm[i][j]*perm[i+1][j]/(perm[i][j]+perm[i+1][j]);
        Kd[i][j]=2*perm[i][j]*perm[i-1][j]/(perm[i][j]+perm[i-1][j]);
        Kl[i][j]=2*perm[i][j]*perm[i][j-1]/(perm[i][j]+perm[i][j-1]);
        Kr[i][j]=2*perm[i][j]*perm[i][j+1]/(perm[i][j]+perm[i][j+1]);
      }
    }
  }

  for (int i=1;i<NY-1;i++)//initializing beta-alpha field
  {
    for (int j=1;j<NX-1;j++)
    {
      if (Kl[i][j] == 0) {
        betal[i][j]=0;
      } else {
        betal[i][j]=cfactor*h/Kl[i][j];
      }
      if (Kr[i][j] == 0) {
        betar[i][j]=0;
      } else {
        betar[i][j]=cfactor*h/Kr[i][j];
      }
      if (Kd[i][j] == 0) {
        betad[i][j]=0;
      } else {
        betad[i][j]=cfactor*h/Kd[i][j];
      }
      if (Ku[i][j] == 0) {
        betau[i][j]=0;
      } else {
        betau[i][j]=cfactor*h/Ku[i][j];
      }
    }
  }
}

/* Prep Data to send */
void ver_prep(double q_data[NY][NX], double l_data[NY][NX], double sendvec[2*NY], int col) {
  for (int i = 0; i<NY; i++) {
    sendvec[i] = q_data[i][col];
  }
  for (int i = 0; i<NY; i++) {
    sendvec[i+NY] = l_data[i][col];
  }
}

void hor_prep(double q_data[NY][NX], double l_data[NY][NX], double sendvec[2*NX], int row) {
  for (int i = 0; i<NX; i++) {
    sendvec[i] = q_data[row][i];
  }
  for (int i = 0; i<NX; i++) {
    sendvec[i+NY] = l_data[row][i];
  }
}

/* Set Buffer Data Functions */
void set_buffer_data_0(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right */
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_1(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right */
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_2(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_3(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right*/
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_4(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right */
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_5(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to bottom */
  MPI_Send(send_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from botton */
  MPI_Recv(receive_down, NX*2, MPI_DOUBLE, block_ID + nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_6(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right */
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_7(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to right */
  MPI_Send(send_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from right */
  MPI_Recv(receive_right, NY*2, MPI_DOUBLE, block_ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void set_buffer_data_8(double send_left[2*NY], double send_right[2*NY],
   double send_down[2*NX], double send_up[2*NX], double receive_left[2*NY],
   double receive_right[2*NY], double receive_down[2*NX],
   double receive_up[2*NX], int block_ID, int nsx) {
  /* exchange info - send to left */
  MPI_Send(send_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD);
  /* exchange info - receive from left */
  MPI_Recv(receive_left, NY*2, MPI_DOUBLE, block_ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  /* exchange info - send to top */
  MPI_Send(send_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD);
  /* exchange info - receive from top */
  MPI_Recv(receive_up, NX*2, MPI_DOUBLE, block_ID - nsx, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Update functions
void ld_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];
  double ql_oldval, ll_oldval, qr_star, lr_star;
  double qd_oldval, ld_oldval, qu_star, lu_star;
  double bdl_type, bdl_val;
  double bdd_type, bdd_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_oldval  = ql_old[i][j];
  ll_oldval = ll_old[i][j];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];
  qd_oldval = qd_old[i][j];
  ld_oldval = ld_old[i][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];

  bdl_type = bd_condl[i][1];
  bdl_val = bd_condl[i][0];
  bdd_type = bd_condd[i][1];
  bdd_val = bd_condd[i][0];

  double chi = (2*k)/h;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Ar = chi*1/(1+beta_r*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1.0 - bdl_type)*bdl_val + (1.0 - bdd_type)*bdd_val;
  double neu = bdl_type*bdl_val + bdd_type*bdd_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Br + Bu)/(Ar + Au + chi*(2.0-bdl_type-bdd_type));

  /* update flux */
  double ql_new = -chi*(bdl_val - p_new)*(1-bdl_type) + bdl_type*bdl_val;
  double qr_new = Ar*p_new - Br;
  double qd_new = -chi*(bdd_val - p_new)*(1-bdd_type) + bdd_type*bdd_val;
  double qu_new = Au*p_new - Bu;

  /* update lagrangians */
  ll[i][j] = bdl_val*(1-bdl_type) + 1/chi*(chi*p_new - bdl_val)*bdl_type;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = bdd_val*(1-bdd_type) + 1/chi*(chi*p_new - bdd_val)*bdd_type;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;
}

void rd_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];

  double ql_star, ll_star, qr_oldval, lr_oldval;
  double qd_oldval, ld_oldval, qu_star, lu_star;
  double bdr_type, bdr_val;
  double bdd_type, bdd_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star  = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_oldval = ql_old[i][j];
  lr_oldval = lr_old[i][j];
  qd_oldval = qd_old[i][j];
  ld_oldval = ld_old[i][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];

  bdr_type = bd_condr[i][1];
  bdr_val = bd_condr[i][0];
  bdd_type = bd_condd[i][1];
  bdd_val = bd_condd[i][0];

  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Al = chi*1/(1+beta_l*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1 - bdr_type)*bdr_val+(1-bdd_type)*(bdd_val);
  double neu = bdr_type*bdr_val + bdd_type*bdd_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Bl + Bu)/(Al + Au + chi*(2-bdr_type-bdd_type));

  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = -chi*(bdr_val - p_new)*(1-bdr_type) + bdr_type*bdr_val;
  double qd_new = -chi*(bdd_val - p_new)*(1-bdd_type) + bdd_type*bdd_val;
  double qu_new = Au*p_new - Bu;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = bdr_val*(1-bdr_type) + 1/chi*(chi*p_new - bdr_val)*bdr_type;
  ld[i][j] = bdd_val*(1-bdd_type) + 1/chi*(chi*p_new - bdd_val)*bdd_type;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;

}

void lu_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];

  double ql_oldval, ll_oldval, qr_star, lr_star;

  double qd_star, ld_star, qu_oldval, lu_oldval;

  double bdl_type, bdl_val;

  double bdu_type, bdu_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_oldval  = ql_old[i][j];
  ll_oldval = ll_old[i][j];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];

  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_oldval = qu_old[i][j];
  lu_oldval = lu_old[i][j];

  bdl_type = bd_condl[i][1];
  bdl_val = bd_condl[i][0];

  bdu_type = bd_condu[i][1];
  bdu_val = bd_condu[i][0];

  double chi = (2*k)/h;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Ar = chi*1/(1+beta_r*chi);
  double Ad = chi*1/(1+beta_d*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1 - bdl_type)*bdl_val+(1-bdu_type)*(bdu_val);
  double neu = bdl_type*bdl_val + bdu_type*bdu_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Br + Bd)/(Ar + Ad + chi*(2-bdl_type-bdu_type));

  /* update flux */
  double ql_new = -chi*(bdl_val - p_new)*(1-bdl_type) + bdl_type*bdl_val;
  double qr_new = Ar*p_new - Br;
  double qd_new = Ad*p_new - Bd;
  double qu_new = -chi*(bdu_val - p_new)*(1-bdu_type) + bdu_type*bdu_val;

  /* update lagrangians */
  ll[i][j] = bdl_val*(1-bdl_type) + 1/chi*(chi*p_new - bdl_val)*bdl_type;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = bdu_val*(1-bdu_type) + 1/chi*(chi*p_new - bdu_val)*bdu_type;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;
}

void ru_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];

  double ql_star, ll_star, qr_oldval, lr_oldval;

  double qd_star, ld_star, qu_oldval, lu_oldval;

  double bdr_type, bdr_val;

  double bdu_type, bdu_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star  = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_oldval = ql_old[i][j];
  lr_oldval = lr_old[i][j];

  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_oldval = qu_old[i][j];
  lu_oldval = lu_old[i][j];

  bdr_type = bd_condr[i][1];
  bdr_val = bd_condr[i][0];
  bdu_type = bd_condu[i][1];
  bdu_val = bd_condu[i][0];

  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Al = chi*1/(1+beta_l*chi);
  double Ad = chi*1/(1+beta_d*chi);
  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1 - bdr_type)*bdr_val+(1-bdu_type)*(bdu_val);
  double neu = bdr_type*bdr_val + bdu_type*bdu_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Bl + Bd)/(Al + Ad + chi*(2-bdr_type-bdu_type));
  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = -chi*(bdr_val - p_new)*(1-bdr_type) + bdr_type*bdr_val;
  double qd_new = Ad*p_new - Bd;
  double qu_new = -chi*(bdu_val - p_new)*(1-bdu_type) + bdu_type*bdu_val;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = bdr_val*(1-bdr_type) + 1/chi*(chi*p_new - bdr_val)*bdr_type;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = bdu_val*(1-bdu_type) + 1/chi*(chi*p_new - bdu_val)*bdu_type;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;
}

void l_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];

  double ql_oldval, ll_oldval, qr_star, lr_star;

  double qd_star, ld_star, qu_star, lu_star;

  double bdl_type, bdl_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_oldval = ql_old[i][j];
  ll_oldval = ll_old[i][j];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];

  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];

  bdl_type = bd_condl[i][1];
  bdl_val = bd_condl[i][0];

  double chi = (2*k)/h;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Ar = chi*1/(1+beta_r*chi);
  double Ad = chi*1/(1+beta_d*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1 - bdl_type)*bdl_val;
  double neu = bdl_type*bdl_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Br + Bd + Bu)/(Ar + Ad + Au + chi*(1-bdl_type));

  /* update flux */
  double ql_new = -chi*(bdl_val - p_new)*(1-bdl_type) + bdl_type*bdl_val;
  double qr_new = Ar*p_new - Br;
  double qd_new = Ad*p_new - Bd;
  double qu_new = Au*p_new - Bu;

  /* update lagrangians */
  ll[i][j] = bdl_val*(1-bdl_type) + 1/chi*(chi*p_new - bdl_val)*bdl_type;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;
}

void r_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];

  double ql_star, ll_star, qr_oldval, lr_oldval;

  double qd_star, ld_star, qu_star, lu_star;

  double bdr_type, bdr_val;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_oldval = qr_old[i][j];
  lr_oldval = lr_old[i][j];

  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];

  bdr_type = bd_condr[i][1];
  bdr_val = bd_condr[i][0];

  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Al = chi*1/(1+beta_l*chi);
  double Ad = chi*1/(1+beta_d*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1 - bdr_type)*bdr_val;
  double neu = bdr_type*bdr_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Bl + Bd + Bu)/(Al + Ad + Au + chi*(1-bdr_type));

  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = -chi*(bdr_val - p_new)*(1-bdr_type) + bdr_type*bdr_val;
  double qd_new = Ad*p_new - Bd;
  double qu_new = Au*p_new - Bu;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = bdr_val*(1-bdr_type) + 1/chi*(chi*p_new - bdr_val)*bdr_type;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;
}

void d_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];
  double ql_star, ll_star, qr_star, lr_star;
  double qd_oldval, ld_oldval, qu_star, lu_star;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];

  qd_oldval = qd_old[i][j];
  ld_oldval = ld_old[i][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];


  double bdd_type, bdd_val;

  bdd_type = bd_condd[i][1];
  bdd_val = bd_condd[i][0];


  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Al = chi*1/(1+beta_l*chi);
  double Ar = chi*1/(1+beta_r*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1-bdd_type)*(bdd_val);
  double neu = bdd_type*bdd_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Bl + Br + Bu)/(Al + Ar + Au + chi*(1-bdd_type));

  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = Ar*p_new - Br;
  double qd_new = -chi*(bdd_val - p_new)*(1-bdd_type) + bdd_type*bdd_val;
  double qu_new = Au*p_new - Bu;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = bdd_val*(1-bdd_type) + 1/chi*(chi*p_new - bdd_val)*bdd_type;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;

}

void u_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];
  double ql_star, ll_star, qr_star, lr_star;
  double qd_star, ld_star, qu_oldval, lu_oldval;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];

  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_oldval = qu_old[i][j];
  lu_oldval = lu_old[i][j];


  double bdu_type, bdu_val;

  bdu_type = bd_condu[i][1];
  bdu_val = bd_condu[i][0];

  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Al = chi*1/(1+beta_l*chi);
  double Ar = chi*1/(1+beta_r*chi);
  double Ad = chi*1/(1+beta_d*chi);

  /* if Dirichlet, then type is 0, else type is 1 */

  double dir = (1-bdu_type)*(bdu_val);
  double neu = bdu_type*bdu_val;

  /* Solve */

  /* pressure */
  double p_new = (fi*h - neu + chi*dir + Bl + Br + Bd)/(Al + Ar + Ad + chi*(1-bdu_type));

  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = Ar*p_new - Br;
  double qd_new = Ad*p_new - Bd;
  double qu_new = -chi*(bdu_val - p_new)*(1-bdu_type) + bdu_type*bdu_val;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = bdu_val*(1-bdu_type) + 1/chi*(chi*p_new - bdu_val)*bdu_type;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;

}

void c_update(int i, int j) {
  /* local variables */
  double beta_l, beta_r, beta_d, beta_u;
  double fi = f[i][j];
  double k = perm[i][j];
  double ql_star, ll_star, qr_star, lr_star;
  double qd_star, ld_star, qu_star, lu_star;

  beta_l = betal[i][j];
  beta_r = betar[i][j];
  beta_d = betad[i][j];
  beta_u = betau[i][j];

  ql_star  = qr_old[i][j-1];
  ll_star = lr_old[i][j-1];
  qr_star = ql_old[i][j+1];
  lr_star = ll_old[i][j+1];
  qd_star = qu_old[i-1][j];
  ld_star = lu_old[i-1][j];
  qu_star = qd_old[i+1][j];
  lu_star = ld_old[i+1][j];

  double chi = (2*k)/h;
  double Bl = chi*1/(1+beta_l*chi)*beta_l*ql_star + chi*1/(1+beta_l*chi)*ll_star;
  double Br = chi*1/(1+beta_r*chi)*beta_r*qr_star + chi*1/(1+beta_r*chi)*lr_star;
  double Bd = chi*1/(1+beta_d*chi)*beta_d*qd_star + chi*1/(1+beta_d*chi)*ld_star;
  double Bu = chi*1/(1+beta_u*chi)*beta_u*qu_star + chi*1/(1+beta_u*chi)*lu_star;
  double Al = chi*1/(1+beta_l*chi);
  double Ar = chi*1/(1+beta_r*chi);
  double Ad = chi*1/(1+beta_d*chi);
  double Au = chi*1/(1+beta_u*chi);

  /* Solve */

  /* pressure */
  double p_new = (fi*h + Bl + Br + Bd + Bu)/(Al + Ar + Ad + Au);

  /* update flux */
  double ql_new = Al*p_new - Bl;
  double qr_new = Ar*p_new - Br;
  double qd_new = Ad*p_new - Bd;
  double qu_new = Au*p_new - Bu;;

  /* update lagrangians */
  ll[i][j] = beta_l*(ql_new + ql_star) + ll_star;
  lr[i][j] = beta_r*(qr_new + qr_star) + lr_star;
  ld[i][j] = beta_d*(qd_new + qd_star) + ld_star;
  lu[i][j] = beta_u*(qu_new + qu_star) + lu_star;

  ql[i][j] = ql_new; qr[i][j] = qr_new; qd[i][j] = qd_new; qu[i][j] = qu_new;

  p[i][j] = p_new;

}

/* PRINT MATRIX AND OUTPUT FUNCTIONS */

void printmatrix(double a[NY][NX]) {
  for (int i = 0; i < NY;i++){
    for (int j = 0;j<NX;j++) {
      double s = a[i][j];
      if (s<0) {printf("%.3f ", s);} else {
        printf(" %.3f ", s);
      }
    }
    printf("\n");
  }
  printf("\n");
}

int save_to_text(int world_rank,int world_size,int nsx,int nsy,double p[NY][NX]){
  char str[100];
  sprintf(str, "%sfile_%d%s", "results/text_", world_rank,".txt");
  printf("%s\n", str);

  char *filename2 = str;

  // open the file for writing
  FILE *fp2 = fopen(filename2, "w");
  if (fp2 == NULL)
  {
      printf("Error opening the file %s", filename2);
      return -1;
  }
  // write to the text file
  fprintf(fp2, "%d \n", world_size);
  fprintf(fp2, "%d \n", world_rank);
  fprintf(fp2, "%d \n", nsx);
  fprintf(fp2, "%d \n", nsy);
  for (int i = 1; i < NY-1; i++) {
    for (int j = 1; j < NX-1; j++) {
      double s = p[i][j];
      if (s<0) {
        fprintf(fp2, "%f ", s);
      } else {
        fprintf(fp2, " %f ", s);
      }
    }
    fprintf(fp2, "\n");
  }
  // close the file
  fclose(fp2);
  return 0;
}

void writeover(double p_old[NY][NX], double p[NY][NX],
  double ql_old[NY][NX], double qr_old[NY][NX], double qd_old[NY][NX], double qu_old[NY][NX],
  double ll_old[NY][NX], double lr_old[NY][NX], double ld_old[NY][NX], double lu_old[NY][NX],
  double ql[NY][NX], double qr[NY][NX], double qd[NY][NX], double qu[NY][NX],
  double ll[NY][NX], double lr[NY][NX], double ld[NY][NX], double lu[NY][NX]) {
  for (int i = 0; i < NY; i++) {
    for (int j = 0; j < NX; j++) {
      p_old[i][j] = p[i][j];
      ql_old[i][j] = ql[i][j];
      qr_old[i][j] = qr[i][j];
      qd_old[i][j] = qd[i][j];
      qu_old[i][j] = qu[i][j];
      ll_old[i][j] = ll[i][j];
      lr_old[i][j] = lr[i][j];
      ld_old[i][j] = ld[i][j];
      lu_old[i][j] = lu[i][j];
    }
  }
}

void convergence(int k, double a[NY][NX], double b[NY][NX], int flag, double myerror[2], int world_rank) {
  double sum1 = myerror[0];
  double sum2 = myerror[1];
  for (int i = 0; i<NY; i++) {
    for (int j = 0; j<NX; j++) {
      sum1 = sum1 + (a[i][j] - b[i][j])*(a[i][j] - b[i][j]);
      sum2 = sum2 + a[i][j]*a[i][j];
    }
  }
  double num = sqrt(sum1);
  double denum = sqrt(sum2);
  if (flag == 0) {
    printf("<Process %d> <ITER %d> LOCAL ERROR = %.9f \n", world_rank, k, num/denum);
  }
  myerror[0] = sum1;
  myerror[1] = sum2;
}

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  printf("<PROCESS %d> REPORT FOR DUTY <TOTAL RANKS %d>\n", world_rank, world_size);
  maxiter = 1500;
  tol = 1e-5;
  int block_type[2];
  if (world_size != nsx*nsy) {
    fprintf(stderr, "World size must be %d for %s\n", nsx*nsy, argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  set_subdomain_ID(block_ID, world_rank, nsx);
  set_domain_type(block_ID, block_type, nsx, nsy);

  void (*set_buffer_data)(double[2*NY], double[2*NY],
   double[2*NX], double[2*NY], double[2*NY],
   double[2*NY], double[2*NX],
   double[2*NX], int , int );
  switch (block_type[1]) {
    case 0:
      set_buffer_data= set_buffer_data_0;
    break;

    case 1:
      set_buffer_data= set_buffer_data_1;
    break;

    case 2:
      set_buffer_data= set_buffer_data_2;
    break;

    case 3:
      set_buffer_data= set_buffer_data_3;
    break;

    case 4:
      set_buffer_data= set_buffer_data_4;
    break;

    case 5:
      set_buffer_data= set_buffer_data_5;
    break;

    case 6:
      set_buffer_data= set_buffer_data_6;
    break;

    case 7:
      set_buffer_data= set_buffer_data_7;
    break;

    case 8:
      set_buffer_data= set_buffer_data_8;
    break;
  }

  // assign update Functions

  void (*update_ld)(int i, int j);

  void (*update_d)(int i, int j);

  void (*update_rd)(int i, int j);

  void (*update_l)(int i, int j);

  void (*update_c)(int i, int j);

  void (*update_r)(int i, int j);

  void (*update_lu)(int i, int j);

  void (*update_u)(int i, int j);

  void (*update_ru)(int i, int j);

  switch (block_type[1]) {
    case 0:
      update_ld = ld_update; /* 1 */
      update_d = d_update;   /* 2 */
      update_rd = d_update;  /* 3 */
      update_l = l_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = l_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = c_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 1:
      update_ld = d_update;  /* 1 */
      update_d = d_update;   /* 2 */
      update_rd = d_update;  /* 3 */
      update_l = c_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = c_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = c_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 2:
      update_ld = d_update;  /* 1 */
      update_d = d_update;   /* 2 */
      update_rd = rd_update; /* 3 */
      update_l = c_update;   /* 4 */
      update_r = r_update;   /* 5 */
      update_lu = c_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = r_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 3:
      update_ld = l_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = c_update; /* 3 */
      update_l = l_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = l_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = c_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 4:
      update_ld = c_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = c_update; /* 3 */
      update_l = c_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = c_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = c_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 5:
      update_ld = c_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = r_update; /* 3 */
      update_l = c_update;   /* 4 */
      update_r = r_update;   /* 5 */
      update_lu = c_update;  /* 6 */
      update_u = c_update;   /* 7 */
      update_ru = r_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 6:
      update_ld = l_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = c_update; /* 3 */
      update_l = l_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = lu_update;  /* 6 */
      update_u = u_update;   /* 7 */
      update_ru = u_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 7:
      update_ld = c_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = c_update; /* 3 */
      update_l = c_update;   /* 4 */
      update_r = c_update;   /* 5 */
      update_lu = u_update;  /* 6 */
      update_u = u_update;   /* 7 */
      update_ru = u_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;

    case 8:
      update_ld = c_update;  /* 1 */
      update_d = c_update;   /* 2 */
      update_rd = r_update; /* 3 */
      update_l = c_update;   /* 4 */
      update_r = r_update;   /* 5 */
      update_lu = u_update;  /* 6 */
      update_u = u_update;   /* 7 */
      update_ru = ru_update;  /* 8 */
      update_c = c_update;   /* 9 */
    break;
  }

  /* Zero out */
  for (int i = 0; i<NY;i++) {
    for (int j = 0; j < NX; j++) {
      perm[i][j] = 0.0;
      f[i][j] = 0.0;
      p[i][j] = 0.0;
      p_old[i][j] = 0.0;
      bd_condl[i][0] = 0.0;
      bd_condr[i][0] = 0.0;
      bd_condd[i][0] = 0.0;
      bd_condu[i][0] = 0.0;
      bd_condl[i][1] = 0.0;
      bd_condr[i][1] = 0.0;
      bd_condd[i][1] = 0.0;
      bd_condu[i][1] = 0.0;

      ql[i][j] = 0;
      qr[i][j] = 0;
      qd[i][j] = 0;
      qu[i][j] = 0;

      ll[i][j] = 0;
      lr[i][j] = 0;
      ld[i][j] = 0;
      lu[i][j] = 0;

      ql_old[i][j] = 0;
      qr_old[i][j] = 0;
      qd_old[i][j] = 0;
      qu_old[i][j] = 0;

      ll_old[i][j] = 0;
      lr_old[i][j] = 0;
      ld_old[i][j] = 0;
      lu_old[i][j] = 0;

      Kl[i][j] = 0;
      Kr[i][j] = 0;
      Kd[i][j] = 0;
      Ku[i][j] = 0;

      betal[i][j] = 0;
      betar[i][j] = 0;
      betad[i][j] = 0;
      betau[i][j] = 0;
    }
  }

  /* Set the Permiability */
  for (int i = 1; i < NY-1; i++ ) {
    for (int j = 1; j < NX-1; j++ ) {
      perm[i][j] = 1.0;
    }
  }

  switch (block_type[1]) {
    case 0: //note that we must have nonzero border on right and up
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0; //right edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 1://note that we must have nonzero border on left and right and up
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0;//right edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 2://note that we must have nonzero border on left and up
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 3://note that we must have nonzero border on left and bottom and up
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0;//right edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottom edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 4:
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0;//right edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottome edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 5:
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottome edge
      }
      for (int j = 1; j < NX-1; j++ ) {
        perm[NY-1][j] = 1.0;//top edge
      }
    break;
    case 6:
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0;//right edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottome edge
      }
    break;
    case 7:
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int i = 1; i < NY-1; i++ ){
        perm[i][NX-1] = 1.0;//right edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottome edge
      }
    break;
    case 8:
      for (int i = 1; i < NY-1; i++ ){
        perm[i][0] = 1.0;//left edge
      }
      for (int j = 1; j < NX-1; j++ ){
        perm[0][j] = 1.0;//bottome edge
      }
    break;
  }

  /* Set up the betas */
  cfactor = 13.0;
  h = 1.0/(((double)NX-2)*nsx);

  set_betas( Kl, Kr, Kd, Ku, perm, betal, betar, betad, betau, cfactor, h);

  /* Set the boundary conditions */
  switch (block_type[1]) {
    case 0: //note that we must have nonzero border on right and up
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condl[i][0] = 1.0;
        bd_condl[i][1] = 0.0;//Dirichlet is zero
      }
      for (int j = 1; j < NX-1 ; j++ ) {
        bd_condd[j][0] = 0.0;
        bd_condd[j][1] = 1.0;
      }
    break;
    case 1://note that we must have nonzero border on left and right and up
      for (int j = 1; j < NX-1 ; j ++ ) {
        bd_condd[j][0] = 0.0;
        bd_condd[j][1] = 1.0;
      }
    break;
    case 2://note that we must have nonzero border on left and up
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condr[i][0] = -1.0;
        bd_condr[i][1] = 0.0;
      }
      for (int j = 1; j < NX-1 ; j ++ ) {
        bd_condd[j][0] = 0.0;
        bd_condd[j][1] = 1.0;//Dirichlet is zero
      }
    break;
    case 3://note that we must have nonzero border on left and bottom and up
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condl[i][0] = 1.0;
        bd_condl[i][1] = 0.0;//Dirichlet is zero
      }
    break;
    case 4:
    break;
    case 5:
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condr[i][0] = -1.0;
        bd_condr[i][1] = 0.0;//Dirichlet is zero
      }
    break;
    case 6:
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condl[i][0] = 1.0;
        bd_condl[i][1] = 0.0;//Dirichlet is zero
      }
      for (int j = 1; j < NX-1 ; j++ ) {
        bd_condu[j][0] = 0.0;
        bd_condu[j][1] = 1.0;//Dirichlet is zero
      }
    break;
    case 7:
      for (int j = 1; j < NX-1 ; j++ ) {
        bd_condu[j][0] = 0.0;
        bd_condu[j][1] = 1.0;//Dirichlet is zero
      }
    break;
    case 8:
      for (int i = 1; i < NY-1 ; i++ ) {
        bd_condr[i][0] = -1.0;
        bd_condr[i][1] = 0.0;//Dirichlet is zero
      }
      for (int j = 1; j < NX-1 ; j++ ) {
        bd_condu[j][0] = 0.0;
        bd_condu[j][1] = 1.0;//Dirichlet is zero
      }
    break;
  }
  printf("<PROCESS %d> block_ID [%d,%d], BLOCK_TYPE %d:>\n", world_rank, block_ID[0],block_ID[1],block_type[1]);
  // enter the loop
  if (world_rank == 8) {
    printf("<PROCESS %d> block_ID [%d,%d], BLOCK_TYPE %d:>\n", world_rank, block_ID[0],block_ID[1],block_type[1]);
    printf("<PROCESS %d> Permiability:\n", world_rank);
    printmatrix(perm);
    printf("<PROCESS %d> BetaL:\n", world_rank);
    printmatrix(betal);
    printf("<PROCESS %d> BetaR:\n", world_rank);
    printmatrix(betar);
    printf("<PROCESS %d> BetaD:\n", world_rank);
    printmatrix(betad);
    printf("<PROCESS %d> BetaU:\n", world_rank);
    printmatrix(betau);
    printf("<PROCESS %d> bd_l val:\n", world_rank);
    for (int i = 0; i<NY; i++) {
      double s = bd_condl[i][0];
      if (s<0) {
        printf("%f ", s);
      } else {
        printf(" %f ", s);
      }
      printf("\n");
    }
    printf("<PROCESS %d> bd_l type:\n", world_rank);
    for (int i = 0; i<NY; i++) {
      double s = bd_condl[i][1];
      if (s<0) {
        printf("%f ", s);
      } else {
        printf(" %f ", s);
      }
      printf("\n");
    }
    printf("<PROCESS %d> bd_U:\n", world_rank);
    for (int j = 0; j<NX; j++) {
      double s = bd_condu[j][0];
      if (s<0) {
        printf("%f ", s);
      } else {
        printf(" %f ", s);
      }
    }
    printf("\n");
  }

  for (int k = 0; k<maxiter; k++) {
    //grid update
    int i = 1;
    int j = 1;
    (*update_ld)(i, j);

    i = 1;
    j = NX-2;
    (*update_rd)(i, j);

    i = NY-2;
    j = 1;
    (*update_lu)(i, j);

    i = NY-2;
    j = NX-2;
    (*update_ru)(i, j);

    i = 1;
    for (j = 2; j<NX-2; j++) {
      (*update_d)(i, j);
    }

    i = NY-2;
    for (j = 2; j<NX-2; j++) {
      (*update_u)(i, j);
    }

    j = 1;
    for (i = 2; i<NY-2; i++) {
      (*update_l)(i, j);
    }

    j = NX-2;
    for (i = 2; i<NY-2; i++) {
      (*update_r)(i, j);
    }

    for (i = 2; i<NY-2; i++) {
      for (j = 2; j<NX-2; j++) {
        (*update_c)(i, j);
      }
    }

    // Check convergence
    int flag = 1;//if flag is zero print local error
    double myerror[2] = {0,0};
    convergence(k, p, p_old, flag, myerror, world_rank);
    double globalerror_num[R], globalerror_denum[R];//is number of processes
    MPI_Allgather(&myerror[0], 1, MPI_DOUBLE, globalerror_num, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    MPI_Allgather(&myerror[1], 1, MPI_DOUBLE, globalerror_denum, 1, MPI_DOUBLE,
            MPI_COMM_WORLD);
    double globalerror;
    double num = 0;
    double denum = 0;
    for (int y = 0; y<R; y++) {
      num = num + globalerror_num[y];
      denum = denum + globalerror_denum[y];
    }
    globalerror = sqrt(num/denum);

    if (globalerror < tol) {
      if (world_rank == 0) {
        printf("<PROCESS %d>: CONVERGED [tol = %f] with iter = %d \n", world_rank, tol, k);
      }
      break;
    }

    if (world_rank == 0) {
      printf("<PROCESS %d> ITER %d: globalerror %f:\n", world_rank, k,globalerror);
      //printmatrix(p);

    }

    // Send and receive
    double send_left[2*NY], send_right[2*NY];
    double send_down[2*NX], send_up[2*NX];
    double receive_left[2*NY], receive_right[2*NY];
    double receive_down[2*NX], receive_up[2*NX];

    int col_left = 1;
    int col_right = NX-2;
    int row_down = NY-2;
    int row_up = 1;

    ver_prep(ql,ll,send_left,col_left);
    ver_prep(qr,lr,send_right,col_right);
    hor_prep(qd,ld,send_up,row_up);
    hor_prep(qu,lu,send_down,row_down);

    (*set_buffer_data)(send_left,send_right,send_down,send_up,
      receive_left, receive_right, receive_down, receive_up, world_rank, nsx);

    for (int i = 0; i<NY; i++) {
      qr[i][0] = receive_left[i];
      lr[i][0] = receive_left[NY+i];
      ql[i][NX-1] = receive_right[i];
      ll[i][NX-1] = receive_right[NY+i];
    }

    for (int j = 0; j<NX; j++) {
      qd[NY-1][j] = receive_down[j];
      ld[NY-1][j] = receive_down[NX+j];
      qu[0][j] = receive_up[j];
      lu[0][j] = receive_up[NX+j];
    }

    //reset
    writeover(p_old, p,
    ql_old, qr_old, qd_old, qu_old,
    ll_old, lr_old, ld_old, lu_old,
    ql, qr, qd, qu,
    ll, lr, ld, lu);

  }
  // Print a matrix out to respective file for all processes
  int x = save_to_text(world_rank,world_size,nsx,nsy,p);


  //check convergence

  MPI_Finalize();
  return 0 ;
}
