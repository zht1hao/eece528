all: test_gen baseline cudaV Pthreads
	@echo ===================================================================
	@echo  This is Haotians part of the group work
	@echo	including a base case, a cuda version and a test case generator.
	@echo	to use the baseline, cuda version and the gernerator, 
	@echo  please do the ./test_gen first to generate some test cases.
	@echo	           ./baseline_gcc bdd_test1.txt bdd_test2.txt
	@echo		  ./cudaV bdd_test1.txt bdd_test.txt
	@echo  You can find the Yu Lei part under his folder.
	@echo ===================================================================

test_gen:
	gcc ./testcase_gen.c -o test_gen
baseline:
	gcc ./Haotian_base.c -lm -O3 -o ./baseline_gcc
	nvcc ./Haotian_base.cu -o ./baseline_cuda
cudaV:
	nvcc ./Haotian_cuda.cu -O3  -o ./cudaV
Pthreads:
	gcc project_pthread.c -lpthread -O2 -o GEP	
clean:
	rm -f test_gen
	rm -f bdd_test.txt
	rm -f bdd_test2.txt
	rm -f baseline_gcc
	rm -f baseline_cuda
	rm -f cudaV

.phony:clean
