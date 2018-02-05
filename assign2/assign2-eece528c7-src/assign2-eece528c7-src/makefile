GNUCompilerCollection = gcc -pthread
Flags = -Wall -Wextra -Ofast

all: Compilation_custom

clean:
        rm -rf *.o *.out

.PHONY: all clean clean_all

Compilation_o3:
        $(GNUCompilerCollection) $(Flags) -O3 -o gaussianelimination_sequentialo3.out gaussianelimination_sequential.c
        $(GNUCompilerCollection) $(Flags) -O3 -o gaussianelimination_pthreadso3.out gaussianelimination_pthreads.c

CUSTOM_FLAGS =  -Wall -Wpedantic -Wextra -O3 -ftree-vectorize -m64 -march=native -std=c99 -pthread
EXTRA_FLAGS = -funroll-loops

Compilation_custom:

        $(GNUCompilerCollection) $(Flags) $(CUSTOM_FLAGS) $(EXTRA_FLAGS) -o gaussianelimination_pthreads.out gaussianelimination_pthreads.c
        $(GNUCompilerCollection) $(Flags) $(CUSTOM_FLAGS) $(EXTRA_FLAGS) -o gaussianelimination_sequential.out gaussianelimination_sequential.c


.PHONY: Compilation_none Compilation_o1 Compilation_o2 Compilation_o3 Compilation_custom
