trapez_integration: trapez_integration.c
	gcc -o trapez_integration -lm trapez_integration.c

trapez_integration_mpi: trapez_integration_mpi.c
	mpicc -o trapez_integration_mpi -lm trapez_integration_mpi.c
trapez_integration_mpi_solved: trapez_integration_mpi_solved.c
	mpicc -o trapez_integration_mpi_solved -lm trapez_integration_mpi_solved.c

clean:
	rm trapez_integration trapez_integration_mpi
