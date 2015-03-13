
/*
 *Library for operation with Ordered Binary Decision Diagrams
 *
 *Author: Volodymyr Mihav
 *The Kirovohrad Volodymyr Vynnychenko State Pedagogical University
 *Date: Februar 2015
 */

#ifndef _BDDCORE_H_
#define _BDDCORE_H_

#ifdef __cplusplus
extern "C"{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if defined (_WIN32) || defined (WIN32)
	#define __windows__
	#include <windows.h>
#elif defined (linux) || defined (__linux)
	#define __linux__
	#include <pthread.h>
#endif

#define LOG_MEMPAGE_SIZE			10	// 0 < LOG_MEMPAGE_SIZE < LOG_MEMCONTAINER_MAXSIZE
#define LOG_MEMCONTAINER_MAXSIZE	31	// 1<<LOG_MEMCONTAINER_MAXSIZE - maximal value signed lo, hi, aux field
	// possible values:	31 - sizeof(lo)=4 (long)

#define MEMPAGE_SIZE		(1<<LOG_MEMPAGE_SIZE)
#define MAX_PAGES			(1<<(LOG_MEMCONTAINER_MAXSIZE-LOG_MEMPAGE_SIZE))

#define node_(x, bdd)			(bdd->container.page[(x)>>LOG_MEMPAGE_SIZE]->node[(x)&(MEMPAGE_SIZE-1)])
#define template_(x, container)		((container).page[(x)>>LOG_MEMPAGE_SIZE]->temp[(x)&(MEMPAGE_SIZE-1)])
#define treeNode_(x)				(tree->node[(x)])

#define varletter_(x)			((x)&1?'A'+(((x)>>1)&31):'a'+(((x)>>1)&31)) //31 - (1<<5)-1
#define varindex_(x)			(x>>6)
#define varhash_(letter, index)	((unsigned short)(((index)<<6)+(letter<'a'?1+((letter-'A')<<1):0+((letter-'a')<<1))))

#pragma pack(push, 2) 

	typedef struct BddNodeType{
		long lo;
		long hi;
		long aux;
		unsigned short v;
		unsigned short ref;
	} BddNode;

	typedef struct PageType{
		long availableNode;
		long availableNodeCount;
		BddNode node[MEMPAGE_SIZE];
	} Page;

	typedef struct MemoryContainerType{
		Page **page;
		long pageCount;
	}MemoryContainer;

	typedef struct BddType{
		unsigned long root; // номер вузла
		unsigned long nodesCount;
		unsigned long terminals[2];
		unsigned short * orderTable;
		unsigned short orderTableSize;
		MemoryContainer container;
	} Bdd;

	typedef struct TemplateType{
		long l;
		long h;
		long left;
		long right;
	}Template;

	typedef struct TemplatePageType{
		long availableTemplateCount;
		Template temp[MEMPAGE_SIZE];
	} TemplatePage;

	typedef struct TemplateContainerType{
		TemplatePage **page;
		long pageCount;
	}TemplateContainer;

#ifdef __windows__
	typedef struct BoolTreeNodeType{
		volatile Bdd * volatile bdd;
		HANDLE thread;
		volatile long parent;
		long left;
		long right;
		unsigned short index;
		unsigned char op;
	}BTreeNode;

	typedef struct BoolTreeContainerType{
		BTreeNode *node;
		CRITICAL_SECTION critSec;
		volatile long availableNode;
	}BTreeContainer;
#endif
#ifdef __linux__
	typedef struct BoolTreeNodeType{
		volatile Bdd * volatile bdd;
		pthread_t thread;
		volatile long parent;
		long left;
		long right;
		unsigned short index;
		unsigned char op;
	}BTreeNode;

	typedef struct BoolTreeContainerType{
		BTreeNode *node;
		pthread_mutex_t mutex;
		volatile long availableNode;
	}BTreeContainer;
#endif

	typedef struct BoolTreeSubContainerType{
		BTreeContainer *tree;
		long index;
	}BTreeSubContainer;
#pragma pack(pop)


	extern long allocateNode(Bdd *bdd);//create new BDD node
	long createPage(Bdd *bdd);//increases the BDD size (it isn't intended for a call by the user)
	extern void freeNode(long node, Bdd *bdd);//places a node in a stack of the free nodes
	extern void freeBdd(Bdd *bdd);
	//{it is necessary for synthesizeBdd()
	long createTemplatePage(TemplateContainer *container);
	long allocateTemplate(TemplateContainer *container, long index);
	void freeTemplates(TemplateContainer *container);
	//}
	//{it is necessary for createBdd()
	long allocateBTreeNode(BTreeContainer *tree);
	void freeBTreeNode(BTreeContainer *tree, long index);
	//}
	long getHash(long l, long r, long logHashTableSize);
	extern long * initializeTraversal(Bdd *bdd);
	extern void deinitializeTraversal(Bdd *bdd, long * head);
	extern void checkBddOrderTable(Bdd *bdd);//deletes variables which aren't present in the BDD
	extern long reduceBdd(Bdd * bdd);//deletes excess nodes
	extern void inverseBdd(Bdd *bdd);//negation the BDD
	extern long swappingVariables(unsigned short upperV, Bdd *bdd, long * head);//exchange of places of variables upperV and upperV+1
	extern unsigned short siftingVariable(unsigned short k, Bdd *bdd, long * head);//search of the best position for a variable k
	extern long makeBddOrder(Bdd *bdd, unsigned short * order, long orderSize);//changes an order of variables
	extern long compareBdd(Bdd *L, Bdd *R);
	extern void compactBdd(Bdd *bdd);//forces out unengaged nodes on top of busy area (can release some pages)
	extern void fprintBdd(FILE * file, Bdd *bdd, long binaryMode);//BDD output to the console or in the file
	extern Bdd* loadBdd(FILE* file);
	extern long optimizeBdd(Bdd *bdd);//search of the best order of variables
	extern long copyBdd(Bdd *target, Bdd *source);
	long findLevel(long f, long g, char op, Bdd *L, Bdd *R);//it is necessary for synthesizeBdd()
	long makeTemplate(long f, long g, long *tbot, long hbase, long logTableSize, TemplateContainer *container, Bdd* L, Bdd *R);//it is necessary for synthesizeBdd()
	extern Bdd* synthesizeBdd(Bdd *L, Bdd *R, char op, long saveL, long saveR);//merge of two BDDs, using the specified operation
	extern Bdd* createBdd(char *boolFunc);//analysis of analytical representation of function and creation of the appropriate BDD

#ifdef __windows__
	DWORD WINAPI createBddThread(LPVOID lpParam);
#endif
#ifdef __linux__
	void* createBddThread(void *lpParam)
#endif

#ifdef __cplusplus
}
#endif

#endif
