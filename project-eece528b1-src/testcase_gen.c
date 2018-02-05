#include<stdio.h>
#include<time.h>
#include<stdlib.h>
int main(){
	FILE *f = fopen("./bdd_test.txt","w");
	FILE *f2 = fopen("./bdd_test2.txt","w");
	
	if(f==NULL){
		printf("Error open file\n");
		exit(1);
	}
	if(f2 == NULL){
		printf("Error open file \n");
		exit(1);
	}
	int zero = -10, one = -11;
	int three,four;
	fprintf(f,"%d %d\n",1000,500);
	fprintf(f,"%d %d %d %d\n",0,32,-10,-11);
	fprintf(f,"%d %d %d %d\n",1,23,0,-11);
	fprintf(f,"%d %d %d %d\n",2,21,-10,0);
	fprintf(f,"%d %d %d %d\n",3,26,2,1);
	srand(time(NULL));
	for(int i=4; i<1000;i++){
		three = rand()%i;
		four = rand()%i;
		if(i % 15 == 0){
			three = -10;
		}
		if(i % 12 == 0){
			four = -11;	
		}
		if(i % 18 ==0){
			three = -11;
		}
		if(i % 20 == 0){
			four = -10;
		}
		
		fprintf(f,"%d %d %d %d\n" ,i, (499-i/2+1),three,four);
			
	}
	fprintf(f2,"%d %d\n",1280,641);
	fprintf(f2,"%d %d %d %d\n",0,32,-10,-11);
	fprintf(f2,"%d %d %d %d\n",1,23,-10,-11);
	fprintf(f2,"%d %d %d %d\n",2,21,1,-10);
	fprintf(f2,"%d %d %d %d\n",3,26,0,1);

	for(int i=4; i<1280;i++){
		three = rand()%i;
		four = rand()%i;
		if(i % 13 == 0){
			three = -10;
		}
		if(i % 21 == 0){
			four = -11;	
		}
		if(i % 31 ==0){
			three = -11;
		}
		if(i % 22 == 0){
			four = -10;
		}
		
		fprintf(f2,"%d %d %d %d\n" ,i, (640-i/2+1),three,four);
			
	}
	fclose(f);	
	fclose(f2);



	return 0;
}
