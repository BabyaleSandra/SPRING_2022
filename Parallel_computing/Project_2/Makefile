CC =gcc
CFLAGS= -Wall -lm
MPICC=mpicc

barkley_serial: barkley_serial.c
	$(MPICC) $(CFLAGS) -o barkley_serial barkley_serial.c

barkley_mpi: barkley_mpi.c
	$(MPICC) $(CFLAGS) -o barkley_mpi barkley_mpi.c

clean:
	rm -f barkley_serial
	rm -f barkley_mpi
