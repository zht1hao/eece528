/*
@EECE528 Project - BDD Parallelization
@Authors: Yu Lei, Haotian Zhang
@Date: 2017/12/3
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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
	bddNode **uniqueNodes;
} applyManager;

void bddTreeInit(bddTree *bdd) {
	bdd->totalNodeNum = 0;
	bdd->totalLevels = 0;
	bddNode *zero, *one;
	zero = (bddNode*)malloc(sizeof(bddNode));
	one = (bddNode*)malloc(sizeof(bddNode));
	zero->index = INFINITY;
	zero->value = 0;
	zero->lowChild = NULL;
	zero->highChild = NULL;
	one->index = INFINITY;
	one->value = 1;
	one->lowChild = NULL;
	one->highChild = NULL;
	bdd->zeroLeaf = zero;
	bdd->oneLeaf = one;
}

void applyManagerInit(applyManager *appMan, int maxNodes) {
	appMan->maxNodeNum = maxNodes;
	appMan->currentSpaceNum = 0;
	bddNode **nodeArray;
	nodeArray = malloc(maxNodes * sizeof(bddNode*));
	appMan->uniqueNodes = nodeArray;
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
	bddNode **nodeArray;
	
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
	nodeArray = malloc(nodeTotal * sizeof(bddNode*));
	bdd->totalNodeNum = nodeTotal;
	bdd->totalLevels = levelTotal;
	while (fgets(linebuf, MAXLINE, f) != NULL) {
		sscanf(linebuf, "%d %d %d %d", &nodeNum, &nodeIndex, &lowC, &highC);
		
		bddNode *newNode;
		newNode = malloc(sizeof(bddNode));
		newNode->index = nodeIndex;
		newNode->value = -1;
		if (lowC ==  -10) {
			newNode->lowChild = bdd->zeroLeaf;
		} else if (lowC == -11) {
			newNode->lowChild = bdd->oneLeaf;
		} else {
			newNode->lowChild = nodeArray[lowC];
		}
		if (highC ==  -10) {
			newNode->highChild = bdd->zeroLeaf;
		} else if (highC == -11) {
			newNode->highChild = bdd->oneLeaf;
		} else {
			newNode->highChild = nodeArray[highC];
		}
		nodeArray[nodeNum] = newNode;
		bdd->topNode = newNode;
	}
	free(nodeArray);
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

bddNode* addBDDs(bddTree *result, bddNode *node1, bddNode *node2, applyManager *appMan) {
	bddNode *left, *right;
	float newNodeIndex;
	int checkNode = 0;
	//printf("CHECKPOINT 1: %.1f, %.1f\n", node1->index, node2->index);
	if (node1->value == 0 && node2->value == 0) {
		return result->zeroLeaf;
	} else if (node1->value == 0 && node2->value == 1) {
		return result->zeroLeaf;
	} else if (node1->value == 1 && node2->value == 0) {
		return result->zeroLeaf;
	} else if (node1->value == 1 && node2->value == 1) {
		return result->oneLeaf;
	}
	
	if (node1->index == node2->index) {
		left = addBDDs(result, node1->lowChild, node2->lowChild, appMan);
		right = addBDDs(result, node1->highChild, node2->highChild, appMan);
	} else if (node1->index < node2->index) {
		left = addBDDs(result, node1->lowChild, node2, appMan);
		right = addBDDs(result, node1->highChild, node2, appMan);
		newNodeIndex = node1->index;
	} else if (node1->index > node2->index) {
		left = addBDDs(result, node1, node2->lowChild, appMan);
		right = addBDDs(result, node1, node2->highChild, appMan);
		newNodeIndex = node2->index;
	}
	
	if (left == right) {
		return left;
	} else {
		checkNode = 0;
		while (checkNode < appMan->currentSpaceNum) {
			if (appMan->uniqueNodes[checkNode]->index == newNodeIndex) {
				if (appMan->uniqueNodes[checkNode]->lowChild == left && appMan->uniqueNodes[checkNode]->lowChild == right) {
					return appMan->uniqueNodes[checkNode];
				}
			}
			//printf("appMan %d %d %d\n", checkNode, appMan->currentSpaceNum, appMan->maxNodeNum);
			checkNode++;
		}
		bddNode *newNode;
		newNode = (bddNode*)malloc(sizeof(bddNode));
		newNode->index = newNodeIndex;
		newNode->value = -1;
		newNode->lowChild = left;
		newNode->highChild = right;
		appMan->uniqueNodes[appMan->currentSpaceNum] = newNode;
		appMan->currentSpaceNum++;
		result->totalNodeNum++;
		//printf("NEW NODE CREATED %.1f %.1f %.1f\n", newNodeIndex, left->index, right->index);
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
	
	bddResult = (bddTree*)malloc(sizeof(bddTree));
	bddTreeInit(bddResult);
	
	applyManager *appMan;
	appMan = (applyManager*)malloc(sizeof(applyManager));
	applyManagerInit(appMan, (int)pow(2, (bdd1->totalLevels + bdd2->totalLevels)));
	
	bddResult->topNode = addBDDs(bddResult, bdd1->topNode, bdd2->topNode, appMan);
	
	printBDD(bdd1);
	printBDD(bdd2);
	printBDD(bddResult);
	
	freeBDD(bdd1);
	freeBDD(bdd2);
	freeBDD(bddResult);
	
	free(appMan->uniqueNodes);
	free(appMan);
	
	return 0;
}