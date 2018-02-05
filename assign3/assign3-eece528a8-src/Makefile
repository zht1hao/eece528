target:
	/usr/local/cuda/bin/nvcc hw3.cu -o "hw3" -DDATA_TYPE=float -DDIM=4096 -O3 
	/usr/local/cuda/bin/nvcc hw3_base.cu -o "hw3_base" -DDATA_TYPE=float -DDIM=4096 -DNTHREADS=1024 -O3 
	/usr/local/cuda/bin/nvcc hw3_cublas.cu -o "hw3_cublas" -DDATA_TYPE=float -DDIM=4096 -lcublas -O3 

