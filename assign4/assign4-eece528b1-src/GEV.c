#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"vbx.h"
#include"vbx_port.h"

#define m_type int

int main(int argc, char** argv){
	#if VBX_SIMULATOR == 1
		vbxsim_init(4,64,256,16,15,4);
	#endif
	int size,i,j,k;
	struct timeval begin2;
	struct timeval end2;
	vbx_timestamp_t begin,end;
	double total_time, total_time2;
	if (argc !=2){
		printf("enter a number to specify the matrix size pls\n for example: ./GEV.elf 4096");
		exit(0);
	}else{
		sscanf(argv[1],"%d",&size);
	}
	m_type *array,*array2,*result,*sumresult,sum;
	
	array = (m_type*)malloc(size*size*sizeof(m_type));
	array2 = (m_type*)malloc(size*size*sizeof(m_type));
	result = (m_type*)malloc(size*size*sizeof(m_type));
	sumresult = (m_type*)malloc(size*size*sizeof(m_type));
	srand(time(NULL));
	
	for(i = 1; i < size; i++){
		for(j = 1; j < size; j++){
			*(array+i*size+j) = rand()%1000;
			*(array2+i*size+j) = rand()%1000;
		}
	}
	
	vbx_timestamp_start();
	begin = vbx_timestamp();
	gettimeofday(&begin2,NULL);

	for(i = 1; i < size; i++){
		for(j = 1; j < size; j++)
			*(array2+j*size +i) = *(array2+i*size+j);	
	}
	
	

	vbx_word_t *va;
	vbx_word_t *vb;
	vbx_word_t *vc;

	vbx_set_vl(size);

	va = (vbx_word_t*)vbx_sp_malloc(size*sizeof(m_type));
	vb = (vbx_word_t*)vbx_sp_malloc(size*sizeof(m_type));
	vc = (vbx_word_t*)vbx_sp_malloc(size*sizeof(m_type));
	
	for (i = 0; i < size; i++){
		vbx_dma_to_vector(va, array+i*size,size*sizeof(m_type));
		for(j = 0; j < size; j++){
			vbx_dma_to_vector(vb,array2+j*size,size*sizeof(m_type));
			vbx_acc(VVW,VMUL,vc,va,vb);
			vbx_dma_to_host(sumresult,vc,size*sizeof(m_type));
		}
//		for(j = 0; j < size; j++){
//			vbx_dma_to_vector(va,sumresult+size*j,size*sizeof(m_type));
//			vbx_acc(VVW,VADD,vc,va,vc);
//		}
//			vbx_dma_to_host(result+size*i,vc,size*sizeof(m_type));
			vbx_sync();	
			
	}
	gettimeofday(&end2,NULL);
	end = vbx_timestamp();
	vbxsim_print_stats_extended();
	vbx_sp_free();
	unsigned freqc;
	freqc = vbx_timestamp_freq();
	total_time = (double) (end-begin)/freqc*1000;
	total_time2 = (double) (end2.tv_sec - begin2.tv_sec) + (double)(end2.tv_usec-begin2.tv_usec)/1000000000;
	printf("total time by vxm: %f\ntotal time by system: %f\n",total_time,total_time2);	
	return 0;
}
