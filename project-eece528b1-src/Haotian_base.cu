/*
@EECE528 Project - BDD Parallelization
@Authors: Yu Lei, Haotian Zhang
@Date: 2017/12/3
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAXNODENUM 	160000
#define MAXLINE 256	/* Maximum length of each input line read. */

typedef struct bddNode_ {
	float index;
	int value;
	struct bddNode_ *lowChild;
	struct bddNode_ *highChild;
} bddNode;

typedef struct bddTree_ {
	int totalNodeNum;
	int totalLevels;
	bddNode *topNode;
	bddNode *zeroLeaf;
	bddNode *oneLeaf;
} bddTree;

typedef struct applyManager_ {
	int maxNodeNum;
	int currentSpaceNum;
} applyManager;
typedef struct pattern_{
	int size;
	float index[MAXNODENUM];
	bddNode* left[MAXNODENUM];
	bddNode* right[MAXNODENUM];
	
}pattern;

pattern patterns;

void bddTreeInit(bddTree *bdd) {
	bddNode *zero,*one;
	bdd->totalNodeNum = 0;
	bdd->totalLevels = 0;
	zero = (bddNode*)malloc(sizeof(bddNode));
	one = (bddNode*)malloc(sizeof(bddNode));
	one->index = INFINITY;
	zero->index = INFINITY;
	zero->value = 0;
	zero->lowChild = NULL;
	zero->highChild = NULL;
	one->value = 1;
	one->lowChild = NULL;
	one->highChild = NULL;
	bdd->zeroLeaf = zero;
	bdd->oneLeaf = one;
}


void applyManagerInit(applyManager *appMan, int maxNodes){
	appMan->maxNodeNum = maxNodes;
	appMan->currentSpaceNum = 0;
}

bddTree* readBDD(char *filename) {
	FILE *f;
	bddTree *bdd;
	int nodeTotal;
	int levelTotal;
	int nodeNum;
	int nodeIndex;
	int lowC;
	int highC;
		
	f = fopen(filename,"r");
	if (!f) {
		fprintf(stderr, "cannot open file \"%s\"\n", filename);
		return NULL;
	}
	
	bdd = (bddTree*)malloc(sizeof(bddTree));
	bddTreeInit(bdd);
	
	char linebuf[MAXLINE];
	
	fgets(linebuf,MAXLINE,f);
	sscanf(linebuf, "%d %d", &nodeTotal, &levelTotal);
	
	bddNode *array[10000];
	
	bdd->totalNodeNum = nodeTotal;
	bdd->totalLevels = levelTotal;
	while (fgets(linebuf, MAXLINE, f) != NULL) {
		sscanf(linebuf, "%d %d %d %d", &nodeNum, &nodeIndex, &lowC, &highC);
		
		bddNode *newNode;
		newNode = (bddNode*)malloc(sizeof(bddNode));
		newNode->index = nodeIndex;
		newNode->value = -1;
		if (lowC ==  -10) {
			newNode->lowChild = bdd->zeroLeaf;
		} else if (lowC == -11) {
			newNode->lowChild = bdd->oneLeaf;
		} else {
			newNode->lowChild = array[lowC];
		}
		if (highC ==  -10) {
			newNode->highChild = bdd->zeroLeaf;
		} else if (highC == -11) {
			newNode->highChild = bdd->oneLeaf;
		} else {
			newNode->highChild = array[highC];
		}
		array[nodeNum] = newNode;
		bdd->topNode = newNode;
	}
	fclose(f);
	return bdd;
}

void printNode(bddNode *node) {
	printf("Node: %f children: \t%f \t%f.\n", node->index, node->lowChild->index, node->highChild->index);
	if (node->lowChild->index != INFINITY) {
		printNode(node->lowChild);
	}
	if (node->highChild->index != INFINITY) {
		printNode(node->highChild);
	}
}

void printBDD(bddTree *bdd) {
	printf("\nPrinting bdd:\n");
	printf("Total nodes in bdd: %d\n", bdd->totalNodeNum);
	
	printNode(bdd->topNode);
}

void recursFree(bddNode *node) {
	if (node->lowChild->index != INFINITY) {
		recursFree(node->lowChild);
	}
	if (node->highChild->index != INFINITY) {
		recursFree(node->highChild);
	}
	free(node);
}

void freeBDD(bddTree *bdd) {
	recursFree(bdd->topNode);
	free(bdd->zeroLeaf);
	free(bdd->oneLeaf);
	free(bdd);
}

// void addNew(int *size) {
	
// }
int check_node(float index,bddNode *left, bddNode* right){
	int i;
	for(i=0;i<patterns.size;i++){
		if(index == patterns.index[i] && left == patterns.left[i] && right == patterns.right[i]){
			return i;	
			
		}
	}
	patterns.index[patterns.size] = index;
	patterns.left[patterns.size] = left;
	patterns.right[patterns.size] = right;
	patterns.size++;
	return 0;	
}
bddNode* applyBDDs(bddTree *result, bddNode *node1, bddNode *node2, applyManager *appMan){
	bddNode *left, *right;
	float newNodeIndex;
	int checkNode = 0;
	
	if(node1->value == 0 && node2->value == 0){
		return result->zeroLeaf;
	}else if(node1->value == 0 && node2->value == 1){
		return result->zeroLeaf;	
	}else if(node1->value == 1 && node2->value == 0){
		return result->zeroLeaf;	
	}else if(node1->value == 1 && node2->value == 1){
		return result->oneLeaf;	
	}
//	printf("node1:%lf node2:%lf",node1->index, node2->index);	
	if(node1->index == node2->index){
		left = applyBDDs(result, node1->lowChild,node2->lowChild,appMan);
		right = applyBDDs(result, node1->highChild,node2->highChild,appMan);
	}else if (node1->index < node2->index){
		left = applyBDDs(result,node1->lowChild,node2,appMan);
		right = applyBDDs(result,node1->highChild,node2,appMan);
		newNodeIndex = node1 -> index;
	}else if (node1->index > node2->index){
		left = applyBDDs(result,node1,node2->lowChild,appMan);
		right = applyBDDs(result,node1,node2->highChild,appMan);
		newNodeIndex == node2 -> index;
	}
//	return result -> oneLeaf;	
	bddNode *newNode;
	newNode = (bddNode*)malloc(sizeof(bddNode));
	if(left == right){
		return left;
	}else{	
		if(checkNode = check_node(newNodeIndex,left,right)){
			newNode->index = patterns.index[checkNode];
			newNode->value = -1;
			newNode->lowChild = patterns.left[checkNode];
			newNode->highChild = patterns.right[checkNode];
		}
		else{
			newNode->index = newNodeIndex;
			newNode->value = -1;
			newNode->lowChild = left;
			newNode->highChild = right;
		}
		return newNode;
	}
	

}

int main(int argc, char* argv[]) {
	bddTree *bdd1, *bdd2;
	bddTree *bddResult;

	
	if (argc !=3) {
		fprintf(stderr,"usage:  a.out file1 file2\n");
		exit(1);
	}
	
    	bdd1 = readBDD(argv[1]);
	bdd2 = readBDD(argv[2]);

//	printBDD(bdd1);
	printf("bdd1 index: %lf and %lf \n",bdd1->topNode->lowChild->index, bdd2->topNode->highChild->index);	
	bddResult = (bddTree*)malloc(sizeof(bddTree));
	bddTreeInit(bddResult);
	
//	printBDD(bdd2);	
	applyManager *appMan;
	appMan = (applyManager*)malloc(sizeof(applyManager));
	applyManagerInit(appMan, (int)pow(2, (bdd1->totalLevels + bdd2->totalLevels)));
	patterns.size = 0;
	
	bddResult->topNode = applyBDDs(bddResult, bdd1->topNode, bdd2->topNode, appMan);
	
//	printBDD(bdd1);
//	printBDD(bdd2);
//	printBDD(bddResult);
	
	free(bdd1);
	free(bdd2);
	free(bddResult);
	
	free(appMan);
	
	return 0;
}
