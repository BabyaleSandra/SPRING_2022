task_2: task_2.c
	gcc -o task_2 -lm task_2.c

std_dev_mpi: std_dev_mpi.c
	mpicc -o std_dev_mpi -lm std_dev_mpi.c

new: task_1.c
	mpicc -o new -lm new.c
clean:
	rm std_dev std_dev_mpi task_1
