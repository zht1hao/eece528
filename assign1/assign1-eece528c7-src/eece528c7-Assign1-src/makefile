GNUCompilerCollection = mpicc
Flags = -Wall -Wextra -Ofast

all: Compilation_custom

clean:
        rm -rf *.o *.out

.PHONY: all clean clean_all

CUSTOM_FLAGS = -O3 -ftree-vectorize -m64 -march=broadwell -mtune=intel
EXTRA_FLAGS = -funroll-loops

Compilation_custom:

        $(GNUCompilerCollection) $(Flags) $(CUSTOM_FLAGS) $(EXTRA_FLAGS) -o gaussianelimination_MPI.out gaussianelimination_MPI.c
        $(GNUCompilerCollection) $(Flags) $(CUSTOM_FLAGS) $(EXTRA_FLAGS) -o gaussianelimination_sequential.out gaussianelimination_sequential.c


.PHONY: Compilation_custom
