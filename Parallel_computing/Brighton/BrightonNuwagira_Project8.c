#include <stdio.h>

#include <mpi.h>

#include <stdlib.h>

#include <math.h>



#define Nx 34

#define Ny 34

#define tolerance 0.000001

#define nsx 2

#define nsy 2

#define N 64





double h= 1.0/(double) N, av,av1,Ar,Al,Au,Ad,Br,Bl,Bu,Bd;

double cfactor = 1.0;

double ql[Ny][Nx],qu[Ny][Nx],qr[Ny][Nx],qd[Ny][Nx];

double ql_old[Ny][Nx],qu_old[Ny][Nx],qr_old[Ny][Nx],qd_old[Ny][Nx];

double ll[Ny][Nx],lu[Ny][Nx],lr[Ny][Nx],ld[Ny][Nx];

double ll_old[Ny][Nx],lu_old[Ny][Nx],lr_old[Ny][Nx],ld_old[Ny][Nx];

double p[Ny][Nx],p_old[Ny][Nx];

double betal[Ny][Nx],betau[Ny][Nx],betar[Ny][Nx],betad[Ny][Nx];

double f[Ny][Nx], perm[Ny][Nx];

double chi[Ny][Nx];

double Ku[Ny][Nx],Kd[Ny][Nx],Kl[Ny][Nx],Kr[Ny][Nx];

double mink=0.0, maxk=0.0,bb,strength=1.0;

double lll=1.0, llr=-1.0;

int world_rank;

int world_size;



int set_subdomain_ID(int rank);

void set_subdomain_type(int rank);

void set_buffer_data_0(int rank);

void set_buffer_data_2(int rank);

void set_buffer_data_6(int rank);

void set_buffer_data_8(int rank);

void interior(int imin,int imax,int jmin,int jmax);

void top_left(int i,int j);

void bottom_left(int i,int j);

void top_right(int i,int j);

void bottom_right(int i,int j );

void top_edge(int jmin,int jmax);

void right_edge(int imin,int imax);

void bottom_edge(int jmin,int jmax);

void left_edge(int imin,int imax);

void update();

void (*set_buffer_data)(int rank);



int block_ID[2];

int block_type;



int main(int argc, char** argv)

{

    MPI_Init(&argc, &argv);



    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);



    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    for(int t=0; t<200;t++)

    {

        if (t ==0)

        {

            FILE* fp = fopen("perm_field_students","r");

                if(fp == NULL)

                {

                    fprintf(stderr,"Unable to open file.Terminating run \n");

                    exit(0);

                }



                for (int j = 0; j<N+2 ;j++)

                {

                    for(int k = 0; k<N+2; k++)

                    {

                        if(k==0 || k==N+1 || j==0 || j==N+1)

                    {

                        perm[k][j]=1.0;

                    }

                    else

                        fscanf(fp,"%lf",&bb);

                        perm[j][k] = exp(strength*bb);

                    }

                }



            //initializing K-alpha field

            set_subdomain_type(world_rank);

            if (block_type == 6)

            {

            for(int i=0;i<Ny;i++)

            {

                for(int j=0;j<Nx;j++)

                {

                    if(i==0 || i==Ny-1 || j==0 || j==Nx-1)

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



            for (int i=1;i<Ny-1;i++)//initializing beta-alpha field

            {

                for (int j=1;j<Nx-1;j++)

                {

                    betau[i][j]=cfactor*h/Ku[i][j];

                    betad[i][j]=cfactor*h/Kd[i][j];

                    betar[i][j]=cfactor*h/Kr[i][j];

                    betal[i][j]=cfactor*h/Kl[i][j];

                }

            }



            }



            if (block_type ==0)

            {



               for(int i=0;i<Ny;i++)

            {

                for(int j=0;j<Nx;j++)

                {

                    if(i==0 || i==Ny-1 || j==0 || j==Nx-1)

                    {

                        Ku[i][j]=0;

                        Kd[i][j]=0;

                        Kl[i][j]=0;

                        Kr[i][j]=0;

                    }

                    else

                    {

                        Ku[i][j]=2*perm[i+Ny-2][j]*perm[i+1+Ny-2][j]/(perm[i+Ny-2][j]+perm[i+1+Ny-2][j]);

                        Kd[i][j]=2*perm[i+Ny-2][j]*perm[i-1+Ny-2][j]/(perm[i+Ny-2][j]+perm[i-1+Ny-2][j]);

                        Kl[i][j]=2*perm[i+Ny-2][j]*perm[i+Ny-2][j-1]/(perm[i+Ny-2][j]+perm[i+Ny-2][j-1]);

                        Kr[i][j]=2*perm[i+Ny-2][j]*perm[i+Ny-2][j+1]/(perm[i+Ny-2][j]+perm[i+Ny-2][j+1]);

                    }

                }

            }



            for (int i=1;i<Ny-1;i++)//initializing beta-alpha field

            {

                for (int j=1;j<Nx-1;j++)

                {

                    betau[i][j]=cfactor*h/Ku[i][j];

                    betad[i][j]=cfactor*h/Kd[i][j];

                    betar[i][j]=cfactor*h/Kr[i][j];

                    betal[i][j]=cfactor*h/Kl[i][j];

                }

            }



            }



            if (block_type == 2)

            {

                    for(int i=0;i<Ny;i++)

            {

                for(int j=0;j<Nx;j++)

                {

                    if(i==0 || i==Ny-1 || j==0 || j==Nx-1)

                    {

                        Ku[i][j]=0;

                        Kd[i][j]=0;

                        Kl[i][j]=0;

                        Kr[i][j]=0;

                    }

                    else

                    {

                        Ku[i][j]=2*perm[i+Ny-2][j+Ny-2]*perm[i+1+Ny-2][j+Ny-2]/(perm[i+Ny-2][j+Ny-2]+perm[i+1+Ny-2][j+Ny-2]);

                        Kd[i][j]=2*perm[i+Ny-2][j+Ny-2]*perm[i-1+Ny-2][j+Ny-2]/(perm[i+Ny-2][j+Ny-2]+perm[i-1+Ny-2][j+Ny-2]);

                        Kl[i][j]=2*perm[i+Ny-2][j+Ny-2]*perm[i+Ny-2][j-1+Ny-2]/(perm[i+Ny-2][j+Ny-2]+perm[i+Ny-2][j-1+Ny-2]);

                        Kr[i][j]=2*perm[i+Ny-2][j+Ny-2]*perm[i+Ny-2][j+1+Ny-2]/(perm[i+Ny-2][j+Ny-2]+perm[i+Ny-2][j+1+Ny-2]);

                    }

                }

            }



            for (int i=1;i<Ny-1;i++)//initializing beta-alpha field

            {

                for (int j=1;j<Nx-1;j++)

                {

                    betau[i][j]=cfactor*h/Ku[i][j];

                    betad[i][j]=cfactor*h/Kd[i][j];

                    betar[i][j]=cfactor*h/Kr[i][j];

                    betal[i][j]=cfactor*h/Kl[i][j];

                }

            }

            }



            if (block_type == 8)

            {

               for(int i=0;i<Ny;i++)

            {

                for(int j=0;j<Nx;j++)

                {

                    if(i==0 || i==Ny-1 || j==0 || j==Nx-1)

                    {

                        Ku[i][j]=0;

                        Kd[i][j]=0;

                        Kl[i][j]=0;

                        Kr[i][j]=0;

                    }

                    else

                    {

                        Ku[i][j]=2*perm[i][j+Ny-2]*perm[i+1][j+Ny-2]/(perm[i][j+Ny-2]+perm[i+1][j+Ny-2]);

                        Kd[i][j]=2*perm[i][j+Ny-2]*perm[i-1][j+Ny-2]/(perm[i][j+Ny-2]+perm[i-1][j+Ny-2]);

                        Kl[i][j]=2*perm[i][j+Ny-2]*perm[i][j-1+Ny-2]/(perm[i][j+Ny-2]+perm[i][j-1+Ny-2]);

                        Kr[i][j]=2*perm[i][j+Ny-2]*perm[i][j+1+Ny-2]/(perm[i][j+Ny-2]+perm[i][j+1+Ny-2]);

                    }

                }

            }



            for (int i=1;i<Ny-1;i++)//initializing beta-alpha field

            {

                for (int j=1;j<Nx-1;j++)

                {

                    betau[i][j]=cfactor*h/Ku[i][j];

                    betad[i][j]=cfactor*h/Kd[i][j];

                    betar[i][j]=cfactor*h/Kr[i][j];

                    betal[i][j]=cfactor*h/Kl[i][j];

                }

            }

            }

            for (int i=0;i<Ny;i++)

            {

                for (int j=0;j<Nx;j++)

                {

                    f[i][j]=0.0;

                    p[i][j]=0.0;

                    qu[i][j]=0.0;

                    qd[i][j]=0.0;

                    qr[i][j]=0.0;

                    ql[i][j]=0.0;

                    lu[i][j]=0.0;

                    ld[i][j]=0.0;

                    lr[i][j]=0.0;

                    ll[i][j]=0.0;

                    chi[i][j]=2*perm[i][j]/h;

                    p_old[i][j]=p[i][j];

                    qu_old[i][j]=qu[i][j];

                    qd_old[i][j]=qd[i][j];

                    qr_old[i][j]=qr[i][j];

                    ql_old[i][j]=ql[i][j];

                    lu_old[i][j]=lu[i][j];

                    ld_old[i][j]=ld[i][j];

                    lr_old[i][j]=lr[i][j];

                    ll_old[i][j]=ll[i][j];

                }

            }

        }

        set_subdomain_type(world_rank);

            switch(block_type)

            {

                case 0:

                    set_buffer_data=&set_buffer_data_0;

                    top_left(Ny-2,1);

                    top_edge(2,Nx-1);

                    left_edge(1,Ny-2);

                    interior(1,Ny-2,2,Nx-1);

                    set_buffer_data(world_rank);

                    update();

                    break;



                case 2:

                    set_buffer_data=&set_buffer_data_2;

                    top_right(Ny-2,Nx-2);

                    top_edge(1,Nx-2);

                    right_edge(1,Ny-2);

                    interior(1,Ny-2,1,Nx-2);

                    set_buffer_data(world_rank);

                    update();

                    break;



                case 6:

                    set_buffer_data=&set_buffer_data_6;

                    bottom_left(1,1);

                    bottom_edge(2,Nx-1);

                    left_edge(2,Ny-1);

                    interior(2,Ny-1,2,Nx-1);

                    set_buffer_data(world_rank);

                    break;



                case 8:

                    set_buffer_data=&set_buffer_data_8;

                    bottom_right(1,Nx-2);

                    bottom_edge(1,Nx-2);

                    right_edge(2,Ny-1);

                    interior(2,Ny-1,1,Nx-2);

                    set_buffer_data(world_rank);

                    break;



            }

            update();



    }

    if (world_rank==0)

            {

                printf("The pressure matrix for case 0 is \n");

                for (int i=Ny-2;i>0;i--)

                    {

                        for (int j=1;j<Nx-1;j++)

                        {

                           printf("%f ",p[i][j]);

                        }

                        printf("\n");

                    }

            }

            if (world_rank==1)

            {

                //update();

                printf("The pressure matrix for case 2 is \n");

                for (int i=i=Ny-2;i>0;i--)

                    {

                        for (int j=1;j<Nx-1;j++)

                        {

                           printf("%f ",p[i][j]);

                        }

                        printf("\n");

                    }

            }

            if (world_rank==2)

            {

                printf("The pressure matrix for case 6 is \n");

                for (int i=Ny-2;i>0;i--)

                    {

                        for (int j=1;j<Nx-1;j++)

                        {

                           printf("%f ",p[i][j]);

                        }

                        printf("\n");

                    }

            }

            if (world_rank==3)

            {

                printf("The pressure matrix for case 8 is \n");

                for (int i=Ny-2;i>0;i--)

                    {

                        for (int j=1;j<Nx-1;j++)

                        {

                           printf("%f ",p[i][j]);

                        }

                        printf("\n");

                    }

            }

    MPI_Finalize();

}



void set_subdomain_type(int rank)

{

    block_ID[0]=rank%nsx;

    block_ID[1]=rank/nsx;

     /* 4 corners */

     if(block_ID[0] == 0 && block_ID[1] == 0)

     block_type = 0;

     else if(block_ID[0] == nsx-1 && block_ID[1] == 0)

     block_type = 2;

     else if(block_ID[0] == nsx-1 && block_ID[1] == nsy-1)

     block_type = 8;

     else if(block_ID[0] == 0 && block_ID[1] == nsy-1)

     block_type = 6;



     /* 4 boundaries */

     else if(block_ID[1] == 0)

     block_type = 1;

     else if(block_ID[0] == nsx-1)

     block_type = 5;

     else if(block_ID[1] == nsy-1)

     block_type = 7;

     else if(block_ID[0] == 0)

     block_type = 3;



     /* interior */

     else

     block_type = 4;

}



void set_buffer_data_0(int rank)

{

    double send_vecr[2*Ny]; //vector for right column values

    for(int i =0;i<Ny;i++)

        {

            send_vecr[i] = qr[i][Nx-2];

            send_vecr[Ny+i] = lr[i][Nx-2];

        }

    MPI_Send(send_vecr,2*Ny,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);



    double rec_vecr[2*Ny]; // vector that recieves right column values

    MPI_Recv(rec_vecr,2*Ny,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        for(int i =0; i<Ny;i++)

        {

            ql[i][Nx-1] = rec_vecr[i];

            ll[i][Nx-1] = rec_vecr [Ny+i];



        }



    double send_vecd[2*Nx]; //vector for bottom row values

    for(int i =0;i<Nx;i++)

        {

            send_vecd[i] = qd[1][i];

            send_vecd[Nx+i] = ld[1][i];

        }

    MPI_Send(send_vecd,2*Nx,MPI_DOUBLE,rank+nsx,0,MPI_COMM_WORLD);



    double rec_vecd[2*Nx]; //vector that recieves bottom row values

    MPI_Recv(rec_vecd,2*Nx,MPI_DOUBLE,rank+nsx,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0;i<Nx;i++)

        {

            qu[0][i] = rec_vecd[i];

            lu[0][i] = rec_vecd[Nx+i];

        }

}

void set_buffer_data_2(int rank)

{

    double send_vecl[2*Ny];//vector for left column values

    for(int i =0;i<Ny;i++)

    {

        send_vecl[i] = ql[i][1];

        send_vecl[Ny+i] = ll[i][1];

    }

    MPI_Send(&send_vecl,2*Ny,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);



    double rec_vecl[2*Ny];// vector that recieves left column values

    MPI_Recv(&rec_vecl,2*Ny,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0; i<Ny;i++)

    {

        qr[i][0] = rec_vecl[i];

        lr[i][0] = rec_vecl[Ny+i];

    }



    double send_vecu[2*Nx];//vector for bottom row values

    for(int i =0;i<Nx;i++)

    {

        send_vecu[i] = qd[1][i];

        send_vecu[Ny+i] = ld[1][i];

    }

    MPI_Send(&send_vecu,2*Nx,MPI_DOUBLE,rank+nsx,0,MPI_COMM_WORLD);



    double rec_vecu[2*Nx];// vector that recieves bottom row values

    MPI_Recv(&rec_vecu,2*Nx,MPI_DOUBLE,rank+nsx,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0; i<Nx;i++)

    {

        qu[0][i] = rec_vecu[i];

        lu[0][i] = rec_vecu[Nx+i];

    }

}



void set_buffer_data_6(int rank)

{

    double send_vecr[2*Ny];//vector for right column values

    for(int i =0;i<Ny;i++)

    {

        send_vecr[i] = qr[i][Nx-2];

        send_vecr[Ny+i] = lr[i][Nx-2];

    }

    MPI_Send(&send_vecr,2*Ny,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);



    double rec_vecr[2*Ny];// vector that recieves left column values

    MPI_Recv(&rec_vecr,2*Ny,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0;i<Ny;i++)

    {

        ql[i][Nx-1] = rec_vecr[i];

        ll[i][Nx-1] = rec_vecr [Ny+i];

    }



    double send_vecu[2*Nx];//vector for top row values

    for(int i =0;i<Nx;i++)

    {

        send_vecu[i] = qu[Ny-2][i];

        send_vecu[Ny+i] = lu[Ny-2][i];

    }

    MPI_Send(&send_vecu,2*Nx,MPI_DOUBLE,rank-nsx,0,MPI_COMM_WORLD);



    double rec_vecd[2*Nx];// vector that top row values

    MPI_Recv(&rec_vecd,2*Nx,MPI_DOUBLE,rank-nsx,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0;i<Ny;i++)

    {

        qd[Ny-1][i] = rec_vecd[i];

        ld[Ny-1][i] = rec_vecd[Ny+i];

    }



}



void set_buffer_data_8(int rank)

{

    double send_vecr[2*Ny];//vector for left column values

    for(int i =0;i<Ny;i++)

    {

        send_vecr[i] = ql[i][1];

        send_vecr[Ny+i] = ll[i][1];

    }

    MPI_Send(&send_vecr,2*Ny,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);



    double rec_vecr[2*Ny];// vector that recieves left column values

    MPI_Recv(&rec_vecr,2*Ny,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0;i<Ny;i++)

    {

        qr[i][0] = rec_vecr[i];

        lr[i][0] = rec_vecr[Ny+i];

    }



    double send_vecu[2*Ny];//vector for top row values

    for(int i =0;i<Nx;i++)

    {

        send_vecu[i] = qu[Ny-2][i];

        send_vecu[Ny+i] = lu[Ny-2][i];

    }

    MPI_Send(&send_vecu,2*Nx,MPI_DOUBLE,rank-nsx,0,MPI_COMM_WORLD);



    double rec_vecu[2*Nx];// vector that recieves top row values

    MPI_Recv(&rec_vecu,2*Nx,MPI_DOUBLE,rank-nsx,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    for(int i =0;i<Nx;i++)

    {

        qd[Ny-1][i] = rec_vecu[i];

        ld[Ny-1][i] = rec_vecu[Ny+i];

    }

}

void interior(int imin,int imax,int jmin,int jmax)

{

    for (int i=imin;i<imax;i++)

    {

        for(int j=jmin;j<jmax;j++)

        {



            Au=chi[i][j]/(1.0+betau[i][j]*chi[i][j]);

            Ad=chi[i][j]/(1.0+betad[i][j]*chi[i][j]);

            Al=chi[i][j]/(1.0+betal[i][j]*chi[i][j]);

            Ar=chi[i][j]/(1.0+betar[i][j]*chi[i][j]);



            Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

            Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

            Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];

            Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];



            p[i][j]=(f[i][j]*h+Bu+Bd+Bl+Br)/(Au+Ad+Al+Ar);



            qr[i][j]=Ar*p[i][j]-Br;

            ql[i][j]=Al*p[i][j]-Bl;

            qu[i][j]=Au*p[i][j]-Bu;

            qd[i][j]=Ad*p[i][j]-Bd;



            lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

            ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

            ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];

            lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];

        }

    }

}



void top_left(int i,int j)

{

    //int i=Ny-2;

    //int j=1;

    Ad=chi[i][j]/(1.0+betad[i][j]*chi[i][j]);

    Ar=chi[i][j]/(1.0+betar[i][j]*chi[i][j]);



    Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

    Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];



    p[i][j]=(f[i][j]*h+chi[i][j]*lll+Bd+Br)/(chi[i][j]+Ad+Ar);



    qr[i][j]=Ar*p[i][j]-Br;

    qd[i][j]=Ad*p[i][j]-Bd;



    ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

    lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];

}



void bottom_left(int i,int j)

{

    //int i=1;

    //int j=1;

    Au=chi[i][j]/(1.0+betau[i][j]*chi[i][j]);

    Ar=chi[i][j]/(1.0+betar[i][j]*chi[i][j]);



    Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

    Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];

    p[i][j]=(f[i][j]*h+chi[i][j]*lll+Bu+Br)/(chi[i][j]+Au+Ar);



    qr[i][j]=Ar*p[i][j]-Br;

    qu[i][j]=Au*p[i][j]-Bu;



    lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

    lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];

}



void top_right(int i,int j)

{

    //i=Ny-2;

    //j=Nx-2;

    Ad=chi[i][j]/(1.0+betad[i][j]*chi[i][j]);

    Al=chi[i][j]/(1.0+betal[i][j]*chi[i][j]);



    Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

    Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];



    p[i][j]=(f[i][j]*h+chi[i][j]*llr+Bl+Bd)/(chi[i][j]+Ad+Al);



    ql[i][j]=Al*p[i][j]-Bl;

    qd[i][j]=Ad*p[i][j]-Bd;



    ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

    ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];



}



void bottom_right(int i,int j)

{

    //int i=1;

    //int j=Nx-2;

    Au=chi[i][j]/(1.0+betau[i][j]*chi[i][j]);

    Al=chi[i][j]/(1.0+betal[i][j]*chi[i][j]);



    Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

    Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];



    p[i][j]=(f[i][j]*h+chi[i][j]*llr+Bu+Bl)/(chi[i][j]+Au+Al);



    ql[i][j]=Al*p[i][j]-Bl;

    qu[i][j]=Au*p[i][j]-Bu;



    lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

    ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];

}



void top_edge(int jmin,int jmax)

{

    //i 1:<Nx-1

    for (int j=jmin;j<jmax;j++)

    {

        int i=Ny-2;

        Ad=chi[i][j]/(1+betad[i][j]*chi[i][j]);

        Al=chi[i][j]/(1+betal[i][j]*chi[i][j]);

        Ar=chi[i][j]/(1+betar[i][j]*chi[i][j]);



        Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

        Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];

        Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];



        p[i][j]=(f[i][j]*h+Bd+Bl+Br)/(Ad+Al+Ar);



        qr[i][j]=Ar*p[i][j]-Br;

        ql[i][j]=Al*p[i][j]-Bl;

        qd[i][j]=Ad*p[i][j]-Bd;



        ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

        ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];

        lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];

    }

}



void right_edge(int imin,int imax)

{

    for (int i=imin;i<imax;i++)

    {

        int j=Nx-2;

        Au=chi[i][j]/(1+betau[i][j]*chi[i][j]);

        Ad=chi[i][j]/(1+betad[i][j]*chi[i][j]);

        Al=chi[i][j]/(1+betal[i][j]*chi[i][j]);



        Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

        Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

        Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];



        p[i][j]=(f[i][j]*h+chi[i][j]*llr+Bu+Bd+Bl)/(chi[i][j]+Au+Ad+Al);



        ql[i][j]=Al*p[i][j]-Bl;

        qu[i][j]=Au*p[i][j]-Bu;

        qd[i][j]=Ad*p[i][j]-Bd;



        lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

        ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

        ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];

    }

}



void bottom_edge(int jmin,int jmax)

{

    for (int j=jmin;j<jmax;j++)

        {

            int i=1;

            Au=chi[i][j]/(1+betau[i][j]*chi[i][j]);

            Al=chi[i][j]/(1+betal[i][j]*chi[i][j]);

            Ar=chi[i][j]/(1+betar[i][j]*chi[i][j]);



            Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

            Bl=Al*betal[i][j]*qr_old[i][j-1]+Al*lr_old[i][j-1];

            Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];



            p[i][j]=(f[i][j]*h+Bu+Bl+Br)/(Au+Al+Ar);

            qr[i][j]=Ar*p[i][j]-Br;

            ql[i][j]=Al*p[i][j]-Bl;

            qu[i][j]=Au*p[i][j]-Bu;



            lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

            ll[i][j]=betal[i][j]*(ql[i][j]+qr_old[i][j-1])+lr_old[i][j-1];

            lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];

        }

}



void left_edge(int imin,int imax)

{

    //i 1:<Ny-2

    for (int i=imin;i<imax;i++)

    {

        int j=1;

        //printf("%f", chi[1][1]);

        Au=chi[i][j]/(1+betau[i][j]*chi[i][j]);

        Ad=chi[i][j]/(1+betad[i][j]*chi[i][j]);

        Ar=chi[i][j]/(1+betar[i][j]*chi[i][j]);



        Bu=Au*betau[i][j]*qd_old[i+1][j]+Au*ld_old[i+1][j];

        Bd=Ad*betad[i][j]*qu_old[i-1][j]+Ad*lu_old[i-1][j];

        Br=Ar*betar[i][j]*ql_old[i][j+1]+Ar*ll_old[i][j+1];



        p[i][j]=(f[i][j]*h+chi[i][j]*lll+Bu+Bd+Br)/(chi[i][j]+Au+Ad+Ar);

        //printf("%f", p[i][j]);



        qr[i][j]=Ar*p[i][j]-Br;

        qu[i][j]=Au*p[i][j]-Bu;

        qd[i][j]=Ad*p[i][j]-Bd;



        lu[i][j]=betau[i][j]*(qu[i][j]+qd_old[i+1][j])+ld_old[i+1][j];

        ld[i][j]=betad[i][j]*(qd[i][j]+qu_old[i-1][j])+lu_old[i-1][j];

        lr[i][j]=betar[i][j]*(qr[i][j]+ql_old[i][j+1])+ll_old[i][j+1];



    }



}

void update()

{

    for (int i=0;i<Ny;i++)

        {

            for (int j=0;j<Nx;j++)

            {

                p_old[i][j]=p[i][j];

                qu_old[i][j]=qu[i][j];

                qd_old[i][j]=qd[i][j];

                qr_old[i][j]=qr[i][j];

                ql_old[i][j]=ql[i][j];

                lu_old[i][j]=lu[i][j];

                ld_old[i][j]=ld[i][j];

                lr_old[i][j]=lr[i][j];

                ll_old[i][j]=ll[i][j];

            }

        }

}
