NVCC	= nvcc
NVCCFLAGS	= 
LIBS	= 
all	: cudaAssign3 cudaAssign3_O cudaAssign3_O1 cudaAssign3_O2 cudaAssign3_O3 cudaAssign3_O4 
clean	:
	rm -f cudaAssign3 cudaAssign3_O cudaAssign3_O1 cudaAssign3_O2 cudaAssign3_O3 cudaAssign3_O4 *.o a.out core

cudaAssign3:
	$(NVCC) -arch=sm_30 -o cudaAssign3 cudaAssign3.cu	$(LIBS)

cudaAssign3_O:
	$(NVCC) -arch=sm_30 -O -o cudaAssign3_O cudaAssign3.cu	$(LIBS)

cudaAssign3_O1:
	$(NVCC) -arch=sm_30 -O1 -o cudaAssign3_O1 cudaAssign3.cu	$(LIBS)

cudaAssign3_O2:
	$(NVCC) -arch=sm_30 -O2 -o cudaAssign3_O2 cudaAssign3.cu      $(LIBS)

cudaAssign3_O3:
	$(NVCC) -arch=sm_30 -O3 -o cudaAssign3_O3 cudaAssign3.cu      $(LIBS)

cudaAssign3_O4:
	$(NVCC) -arch=sm_30 -Xptxas -O3,-v -o cudaAssign3_O4 cudaAssign3.cu      $(LIBS)

