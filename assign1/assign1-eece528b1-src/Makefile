all: baseline Vmpi
	@echo ==========================================================
	@echo to get different size, enter the size after the exe
	@echo for example:
	@echo 		./ge_base 4096
	@echo           mpirun -np32 ge_mpi 4096
	@echo Notice:ge_base is the baseline
	@echo        ge_mpi is the mpi version
	@echo PM: there is a old version of Send/Recv experiment
	@echo         Name: ge_mpi_old
	@echo         Detail:failed experment
	@echo ==========================================================	
baseline: 
	gcc GE.c -O3 -ftree-vectorize -ftree-loop-vectorize -o ge_base
Vmpi:
	mpicc GEM3.c -o ge_mpi -DNDEBUG -O2
Old:
	mpicc GEM.c -o ge_mpi_old -DNDEBUG -O2
