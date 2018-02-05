warmup:warmup.c
	gcc warmup.c -o warmup
warmup_O1:warmup.c
	gcc -O1 warmup.c -o warmup_O1
warmup_O2:warmup.c
	gcc -O2 warmup.c -o warmup_O2
warmup_O3:warmup.c
	gcc -O3 warmup.c -o warmup_O3
warmup_Of:warmup.c
	gcc -O3 -ftree-vectorize -m64 -march=broadwell -mtune=intel -funroll-all-loops -floop-block  warmup.c -o warmup_Of



mpi:mpi.c
	mpicc mpi.c -o mpi
mpi_O1:mpi.c
	mpicc -O1 mpi.c -o mpi_O1
mpi_O2:mpi.c
	mpicc -O2 mpi.c -o mpi_O2
mpi_O3:mpi.c
	mpicc -O3 mpi.c -o mpi_O3
mpi_Of:mpi.c
	mpicc -O3 -ftree-vectorize -m64 -march=broadwell -mtune=intel -funroll-all-loops -floop-block  mpi.c -o mpi_Of