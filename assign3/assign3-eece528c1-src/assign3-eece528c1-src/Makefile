CC      = nvcc
CC2     = gcc
CFLAGS  =
LIBS    =
all     :assign1 assign1_O assign1_O1 assign1_O2 assign1_O3 assign1_O2M \
        assign1_CUDA assign1_CUDA_O assign1_CUDA_O1 assign1_CUDA_O2 assign1_CUDA_O3 assign1_CUDA_O2M

clean:
        rm -f assign1 assign1_O assign1_O1 assign1_O2 assign1_O3 assign1_O2M \
        assign1_CUDA assign1_CUDA_O assign1_CUDA_O1 assign1_CUDA_O2 assign1_CUDA_O3 \
        assign1_CUDA_O2M  *.o core

assign1:
        $(CC2) -o assign1 gaussC.c      $(LIBS)

assign1_O:
        $(CC2) -o assign1_O -O gaussC.c $(LIBS)

assign1_O1:
        $(CC2) -o assign1_O1 -O1 gaussC.c       $(LIBS)

assign1_O2:
        $(CC2) -o assign1_O2 -O2 gaussC.c      $(LIBS)

assign1_O3:
        $(CC2) -o assign1_O3 -O3 gaussC.c      $(LIBS)

assign1_O2M:
        $(CC2) -o assign1_O2M -O2 -ftree-vectorize -m64 -march=broadwell -mtune=intel gaussC.c      $(LIBS)

assign1_CUDA:
        $(CC) -o assign1_CUDA gaussCUDA.cu    $(LIBS)

assign1_CUDA_O:
        $(CC) -O -o assign1_CUDA_O gaussCUDA.cu       $(LIBS)

assign1_CUDA_O1:
        $(CC) -O1 -o assign1_CUDA_O1 gaussCUDA.cu     $(LIBS)

assign1_CUDA_O2:
        $(CC) -O2 -o assign1_CUDA_O2 gaussCUDA.cu      $(LIBS)

assign1_CUDA_O3:
        $(CC) -O3 -o assign1_CUDA_O3 gaussCUDA.cu      $(LIBS)
