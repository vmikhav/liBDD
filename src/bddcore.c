#include "bddcore.h"

long createPage(Bdd *bdd){
	long i, j, offset;
	Page ** tempContainer;

	for (i = 0; i < bdd->container.pageCount; i++){
		if (bdd->container.page[i] == NULL){ break; }
	}
	if (i == bdd->container.pageCount){
		if (bdd->container.pageCount == MAX_PAGES) { return -1; }
		else{
			if (bdd->container.pageCount == 0){
				bdd->container.page = (Page**)malloc(sizeof(Page*));
				bdd->container.pageCount++;
			}
			else{
				tempContainer = (Page**)malloc(sizeof(Page*)*bdd->container.pageCount * 2);
				for (j = 0; j < bdd->container.pageCount; j++){
					tempContainer[j] = bdd->container.page[j];
				}
				bdd->container.pageCount *= 2;
				for (; j < bdd->container.pageCount; j++){
					tempContainer[j] = NULL;
				}
				free(bdd->container.page); bdd->container.page = tempContainer;
			}
		}
	}
	bdd->container.page[i] = (Page*)calloc(1, sizeof(Page));
	bdd->container.page[i]->node[MEMPAGE_SIZE - 1].hi = -1;
	offset = i*MEMPAGE_SIZE;
	for (j = 0; j < MEMPAGE_SIZE - 1; j++){
		bdd->container.page[i]->node[j].hi = offset + j + 1;

	}
	bdd->container.page[i]->availableNode = offset;
	bdd->container.page[i]->availableNodeCount = MEMPAGE_SIZE;
	return i;
}

long allocateNode(Bdd *bdd){
	long i = 0, j = 0;
	for (i = 0; i < bdd->container.pageCount; i++){
		if ((bdd->container.page[i] != NULL) && (bdd->container.page[i]->availableNode != -1)) { break; };
	}
	if (i == bdd->container.pageCount){
		i = createPage(bdd);
		if (i == -1){
			exit(-1073740791);
			// Memory buffer overflow, there are no free pages
		}
	}
	j = bdd->container.page[i]->availableNode;
	bdd->container.page[i]->availableNode = node_(j, bdd).hi;
	bdd->container.page[i]->availableNodeCount--;
	node_(j, bdd).aux = 0;
	node_(j, bdd).v = 1;
	return j;
}

void freeNode(long node, Bdd *bdd){
	long i = 1;
	if (node_(node, bdd).v > 0){
		node_(node, bdd).hi = bdd->container.page[node / MEMPAGE_SIZE]->availableNode;
		node_(node, bdd).v = 0; node_(node, bdd).aux = 0;
		bdd->container.page[node / MEMPAGE_SIZE]->availableNode = node;
		bdd->container.page[node / MEMPAGE_SIZE]->availableNodeCount++;

		if (bdd->container.page[node / MEMPAGE_SIZE]->availableNodeCount == MEMPAGE_SIZE){
			free(bdd->container.page[node / MEMPAGE_SIZE]);
			bdd->container.page[node / MEMPAGE_SIZE] = NULL;
		}
	}
}

void freeBdd(Bdd *bdd){
	long i;
	bdd->nodesCount = 0;
	bdd->orderTableSize = 0;
	bdd->root = 0;
	bdd->terminals[0] = 0; bdd->terminals[1] = 0;

	for (i = 0; i < bdd->container.pageCount; i++){
		if (bdd->container.page[i] != NULL){
			free(bdd->container.page[i]);
		}
	}
	free(bdd->container.page); bdd->container.pageCount = 0;
	free(bdd->orderTable);

	free(bdd);
}

long createTemplatePage(TemplateContainer *container){
	long i, j;
	TemplatePage ** tempContainer;

	for (i = 0; i < container->pageCount; i++){
		if (container->page[i] == NULL){ break; }
	}
	if (i == container->pageCount){
		if (container->pageCount == MAX_PAGES) { return -1; }
		else{
			if (container->pageCount == 0){
				container->page = (TemplatePage**)malloc(sizeof(TemplatePage*));
				container->pageCount++;
			}
			else{
				tempContainer = (TemplatePage**)malloc(sizeof(TemplatePage*)*container->pageCount * 2);
				for (j = 0; j < container->pageCount; j++){
					tempContainer[j] = container->page[j];
				}
				container->pageCount *= 2;
				for (; j < container->pageCount; j++){
					tempContainer[j] = NULL;
				}
				free(container->page); container->page = tempContainer;
			}
		}
	}
	container->page[i] = (TemplatePage*)calloc(1, sizeof(TemplatePage));
	container->page[i]->availableTemplateCount = MEMPAGE_SIZE;
	return i;
}

long allocateTemplate(TemplateContainer *container, long index){
	long i = 0;
	long j = (index >> LOG_MEMPAGE_SIZE) + 1;
	for (i = 0; i < j; i++){
		if ((container->pageCount <= i) || (container->page[i] == NULL)){
			if (createTemplatePage(container) == -1){
				exit(-1073740791);
				// Memory buffer overflow, there are no free pages
			}
		}
	}

	return 1;
}

void freeTemplates(TemplateContainer *container){
	long i;
	for (i = 0; i < container->pageCount; i++){
		if (container->page[i] != NULL){
			free(container->page[i]);
		}
	}
	free(container->page); container->pageCount = 0;
	container->page = NULL;
}

#ifdef __windows__
long allocateBTreeNode(BTreeContainer *tree){
	long i;
	do{
		EnterCriticalSection(&(tree->critSec));
		i = tree->availableNode;
		if (i != -1){ tree->availableNode = treeNode_(tree->availableNode).parent; }
		LeaveCriticalSection(&(tree->critSec));
		if (i == -1){ Sleep(2); }
		else{ break; }
	} while (i == -1);
	return i;
}

void freeBTreeNode(BTreeContainer *tree, long index){
	EnterCriticalSection(&(tree->critSec));
	treeNode_(index).bdd = NULL;
	treeNode_(index).thread = NULL;
	treeNode_(index).parent = tree->availableNode;
	tree->availableNode = index;
	LeaveCriticalSection(&(tree->critSec));
}
#endif
#ifdef __linux__
long allocateBTreeNode(BTreeContainer *tree){
	long i;
	do{
		pthread_mutex_lock(&(tree->mutex));
		i = tree->availableNode;
		if (i != -1){ tree->availableNode = treeNode_(tree->availableNode).parent; }
		pthread_mutex_unlock(&(tree->mutex));
		if (i == -1){ Sleep(2); }
		else{ break; }
	} while (i == -1);
	return i;
}

void freeBTreeNode(BTreeContainer *tree, long index){
	pthread_mutex_lock(&(tree->mutex));
	treeNode_(index).bdd = NULL;
	treeNode_(index).thread = NULL;
	treeNode_(index).parent = tree->availableNode;
	tree->availableNode = index;
	pthread_mutex_unlock(&(tree->mutex));
}
#endif

long getHash(long l, long r, long logHashTableSize){
	//314159257 - A247351 (#9)
	//271828171 - A058814 (#20)
	long i = logHashTableSize / 2;
	return ((314159257 * l + 271828171 * r)&((1 << (16 + i)) - 1)) >> ((16 + i) - logHashTableSize);
}

//prepares the BDD for preorder traversal
//returns the table of pointers on BDD levels
long * initializeTraversal(Bdd *bdd){

	long p, s;
	long * head = (long*)malloc(bdd->orderTableSize*sizeof(long));

	for (p = 0; p<bdd->orderTableSize; p++) { head[p] = -1; }
	s = bdd->root; node_(bdd->root, bdd).aux = ~0; node_(bdd->terminals[0], bdd).aux = ~0; node_(bdd->terminals[1], bdd).aux = ~0;
	while (s != 0){
		p = s;
		s = ~node_(p, bdd).aux;
		node_(p, bdd).aux = head[node_(p, bdd).v];
		head[node_(p, bdd).v] = ~p;
		node_(p, bdd).ref = 0;
		if (node_(node_(p, bdd).lo, bdd).aux >= 0){ node_(node_(p, bdd).lo, bdd).aux = ~s; s = node_(p, bdd).lo; }
		if (node_(node_(p, bdd).hi, bdd).aux >= 0){ node_(node_(p, bdd).hi, bdd).aux = ~s; s = node_(p, bdd).hi; }
	}

	for (s = 0; s < bdd->orderTableSize; s++){
		p = ~head[s];
		while (p){
			node_(node_(p, bdd).lo, bdd).ref++;
			node_(node_(p, bdd).hi, bdd).ref++;
			p = ~node_(p, bdd).aux;
		}
	}
	node_(bdd->root, bdd).ref = 1;

	return head;
}

void deinitializeTraversal(Bdd *bdd, long * head){
	long p, p2, i;
	for (i = 0; i < bdd->orderTableSize; i++){
		p = ~head[i];
		while (p){
			p2 = ~node_(p, bdd).aux;
			node_(p, bdd).aux = 0;
			p = p2;
		}
	}
	free(head);
}

void checkBddOrderTable(Bdd *bdd){
	long i, j, p;
	long *head = initializeTraversal(bdd);
	unsigned short *tempOrderTable;

	for (i = 1, j = 1; i < bdd->orderTableSize; i++, j++){
		if (head[i] == -1){ j--; continue; }
		else if (i != j){
			bdd->orderTable[j] = bdd->orderTable[i];
			p = ~head[i];
			while (p){
				node_(p, bdd).v = (unsigned short)j;
				p = ~node_(p, bdd).aux;
			}
			head[j] = head[i];
		}
	}
	if (i != j){
		if (j > 1){
			tempOrderTable = (unsigned short*)malloc(sizeof(short)* j);
			for (i = 1; i < j; i++){
				tempOrderTable[i] = bdd->orderTable[i];
			}
		}
		else{
			tempOrderTable = NULL; j = 0;
		}
		free(bdd->orderTable); bdd->orderTableSize = j;
		bdd->orderTable = tempOrderTable;
	}

	deinitializeTraversal(bdd, head);
}

//deletes isomorphic subgraphs. returns quantity of the disposed nodes
//returns quantity of the disposed nodes
//not to use between calls initializeTraversal() and deinitializeTraversal()
long reduceBdd(Bdd * bdd){

	//algorithm R (chapter 7.1.4)

	long avail = -1, p = 0, s = 0, v = 0, p2 = 0, q = 0, r = 0;

	long * head = (long*)malloc(bdd->orderTableSize*sizeof(long));
	//R1
	for (p = 0; p<bdd->orderTableSize; p++) { head[p] = -1; }
	s = bdd->root; node_(bdd->root, bdd).aux = ~0; node_(bdd->terminals[0], bdd).aux = ~0; node_(bdd->terminals[1], bdd).aux = ~0;
	while (s != 0){
		p = s;
		s = ~node_(p, bdd).aux;
		node_(p, bdd).aux = head[node_(p, bdd).v];
		head[node_(p, bdd).v] = ~p;
		if (node_(node_(p, bdd).lo, bdd).aux >= 0){ node_(node_(p, bdd).lo, bdd).aux = ~s; s = node_(p, bdd).lo; }
		if (node_(node_(p, bdd).hi, bdd).aux >= 0){ node_(node_(p, bdd).hi, bdd).aux = ~s; s = node_(p, bdd).hi; }
	}
	// R2
	node_(bdd->terminals[0], bdd).aux = 0; node_(bdd->terminals[1], bdd).aux = 0; v = bdd->orderTableSize - 1;
	// R3
	do{
		p = ~head[v]; s = 0;
		while (p){
			p2 = ~node_(p, bdd).aux;
			q = node_(p, bdd).hi; if (node_(q, bdd).lo < 0) { node_(p, bdd).hi = ~node_(q, bdd).lo; }
			q = node_(p, bdd).lo; if (node_(q, bdd).lo < 0) { node_(p, bdd).lo = ~node_(q, bdd).lo; q = node_(p, bdd).lo; }
			if (q == node_(p, bdd).hi){
				node_(p, bdd).lo = ~q;
				node_(p, bdd).hi = avail;
				node_(p, bdd).aux = 0;
				avail = p;
			}
			else if (node_(q, bdd).aux >= 0){
				node_(p, bdd).aux = s;
				s = ~q;
				node_(q, bdd).aux = ~p;
			}
			else {
				node_(p, bdd).aux = node_(~node_(q, bdd).aux, bdd).aux;
				node_(~node_(q, bdd).aux, bdd).aux = p;
			}
			p = p2;
		}
		// R4
		r = ~s; s = 0;
		while (r >= 0){
			q = ~node_(r, bdd).aux; node_(r, bdd).aux = 0;
			if (s == 0) { s = q; }
			else { node_(p, bdd).aux = q; }
			p = q;
			while (node_(p, bdd).aux > 0) p = node_(p, bdd).aux;
			r = ~node_(p, bdd).aux;
		}
		// R5
		p = s; if (p != 0){ q = p; }
		while (p != 0){
			// R6
			s = node_(p, bdd).lo;
			// R7
			do{
				r = node_(q, bdd).hi;
				if (node_(r, bdd).aux >= 0){
					node_(r, bdd).aux = ~q;
				}
				else{
					node_(q, bdd).lo = node_(r, bdd).aux;
					node_(q, bdd).hi = avail; avail = q;
				}
				q = node_(q, bdd).aux;
			} while ((q != 0) && (node_(q, bdd).lo == s));
			// R8
			do{
				if (node_(p, bdd).lo >= 0) { node_(node_(p, bdd).hi, bdd).aux = 0; }
				p = node_(p, bdd).aux;
			} while (p != q);
			// R9
		}
	} while ((v>node_(bdd->root, bdd).v) && (v = v - 1));
	if (node_(bdd->root, bdd).lo < 0) { bdd->root = ~node_(bdd->root, bdd).lo; }

	r = 0;
	while (avail >= 0){
		p = avail;
		avail = node_(avail, bdd).hi;
		freeNode(p, bdd);
		bdd->nodesCount--;
		r++;
	}

	free(head);
	checkBddOrderTable(bdd);
	return r;
}

void inverseBdd(Bdd *bdd){
	bdd->terminals[0] ^= 1;
	node_(0, bdd).lo ^= 1;
	node_(0, bdd).hi ^= 1;

	bdd->terminals[1] ^= 1;
	node_(1, bdd).lo ^= 1;
	node_(1, bdd).hi ^= 1;
}

//swapping of variables with the upperV and upperV+1 indexes
//long * head - result of initializeTraversal()
//returns quantity of the disposed nodes
long swappingVariables(unsigned short upperV, Bdd *bdd, long * head){

	if (upperV == bdd->orderTableSize - 1) return 0;
	if (head == NULL) return 0;

	long a, b, c, d, i, j, k, reserveList = -1, result = 0, uniqueTableSize = 0, logUniqueTableSize = 0;
	long *uniqueTable;

	//adds up quantity of nodes of the upperV level, disunites levels among themselves
	i = ~head[upperV]; j = -1;
	while (i){
		if ((node_(node_(i, bdd).lo, bdd).v != (upperV + 1)) && (node_(node_(i, bdd).hi, bdd).v != (upperV + 1))){
			//the node specifies a sink
			k = node_(i, bdd).aux;
			if (j == -1){ head[upperV] = k; }
			else{ node_(j, bdd).aux = k; }
			node_(i, bdd).aux = reserveList;
			reserveList = ~i;
			node_(i, bdd).v++;
			i = ~k; continue;
		}
		else{
			uniqueTableSize++;
			node_(node_(i, bdd).lo, bdd).ref--;
			node_(node_(i, bdd).hi, bdd).ref--;
		}
		j = i; i = ~node_(i, bdd).aux;
	}

	//creates the table of uniqueness of nodes
	for (i = 2, j = 1; i < uniqueTableSize; i *= 2, j++);
	uniqueTableSize = i;
	uniqueTable = (long*)malloc(sizeof(long)*uniqueTableSize);
	logUniqueTableSize = j;
	for (i = 0; i < uniqueTableSize; i++){ uniqueTable[i] = -1; }

	//adding nodes in the table of uniqueness for a reuse
	i = ~reserveList;
	while (i){
		a = ~node_(i, bdd).aux;
		b = getHash(node_(i, bdd).lo, node_(i, bdd).hi, logUniqueTableSize);
		node_(i, bdd).aux = uniqueTable[b];
		uniqueTable[b] = ~i;
		i = a;
	}

	//checks nodes of the upperV+1 level for possibility of a reuse
	i = ~head[upperV + 1]; reserveList = -1;
	while (i){
		a = ~node_(i, bdd).aux;
		if (node_(i, bdd).ref>0){
			node_(i, bdd).aux = reserveList;
			reserveList = ~i;
		}
		else {
			b = getHash(node_(i, bdd).lo, node_(i, bdd).hi, logUniqueTableSize);
			node_(i, bdd).aux = uniqueTable[b];
			uniqueTable[b] = ~i;
		}
		i = a;
	}
	head[upperV + 1] = -1;

	//swapping of the levels
	i = ~head[upperV]; j = -1;
	while (i){
		//getting ways from a node of the upperV level
		if (node_(node_(i, bdd).lo, bdd).v == (upperV + 1)){
			a = node_(node_(i, bdd).lo, bdd).lo;
			b = node_(node_(i, bdd).lo, bdd).hi;
		}
		else{
			a = node_(i, bdd).lo; b = a;
		}
		if (node_(node_(i, bdd).hi, bdd).v == (upperV + 1)){
			c = node_(node_(i, bdd).hi, bdd).lo;
			d = node_(node_(i, bdd).hi, bdd).hi;
		}
		else{
			c = node_(i, bdd).hi; d = c;
		}


		if (a == c){ node_(i, bdd).lo = a; node_(a, bdd).ref++; }//the lo field of a node specifies a sink
		else{//search of a suitable node for addition of a subfunction
			k = ~uniqueTable[getHash(a, c, logUniqueTableSize)];
			while (k){
				if ((node_(k, bdd).lo == a) && (node_(k, bdd).hi == c)){
					node_(i, bdd).lo = k; node_(k, bdd).ref++;
					k = -1; break;
				}
				k = ~node_(k, bdd).aux;
			}
			if (k > -1){//the suitable node isn't found, creation of the new
				k = allocateNode(bdd); result--;
				node_(k, bdd).v = upperV + 1; node_(k, bdd).ref = 1;
				node_(k, bdd).lo = a; node_(k, bdd).hi = c;
				node_(a, bdd).ref++; node_(c, bdd).ref++;
				node_(k, bdd).aux = uniqueTable[getHash(a, c, logUniqueTableSize)];
				uniqueTable[getHash(a, c, logUniqueTableSize)] = ~k;
				node_(i, bdd).lo = k;
			}
		}

		if (b == d){ node_(i, bdd).hi = b; node_(b, bdd).ref++; }//the lo field of a node specifies a sink
		else{//search of a suitable node for addition of a subfunction
			k = ~uniqueTable[getHash(b, d, logUniqueTableSize)];
			while (k){
				if ((node_(k, bdd).lo == b) && (node_(k, bdd).hi == d)){
					node_(i, bdd).hi = k; node_(k, bdd).ref++;
					k = -1; break;
				}
				k = ~node_(k, bdd).aux;
			}
			if (k > -1){//the suitable node isn't found, creation of the new
				k = allocateNode(bdd); result--;
				node_(k, bdd).v = upperV + 1; node_(k, bdd).ref = 1;
				node_(k, bdd).lo = b; node_(k, bdd).hi = d;
				node_(b, bdd).ref++; node_(d, bdd).ref++;
				node_(k, bdd).aux = uniqueTable[getHash(b, d, logUniqueTableSize)];
				uniqueTable[getHash(b, d, logUniqueTableSize)] = ~k;
				node_(i, bdd).hi = k;
			}
		}
		j = i; i = ~node_(i, bdd).aux;
	}

	//relocation of nodes which can't be changed, to the upperV level
	i = ~reserveList;
	while (i){
		j = ~node_(i, bdd).aux;
		node_(i, bdd).v--;
		node_(i, bdd).aux = head[upperV];
		head[upperV] = ~i;
		i = j;
	}

	//deleting unused nodes
	for (i = 0; i < uniqueTableSize; i++){
		j = ~uniqueTable[i];
		while (j){
			k = ~node_(j, bdd).aux;
			if (node_(j, bdd).ref == 0){
				node_(node_(j, bdd).lo, bdd).ref--;
				node_(node_(j, bdd).hi, bdd).ref--;
				freeNode(j, bdd); result++;
			}
			else{
				node_(j, bdd).aux = head[upperV + 1];
				head[upperV + 1] = ~j;
			}
			j = k;
		}
	}

	//swapping of indexes of variables
	i = bdd->orderTable[upperV];
	bdd->orderTable[upperV] = bdd->orderTable[upperV + 1];
	bdd->orderTable[upperV + 1] = (unsigned short)i;

	free(uniqueTable);
	bdd->nodesCount -= result;

	return result;
}

//search of an optimum location for a variable k
//long * head - result of initializeTraversal()
//returns new location of a variable
unsigned short siftingVariable(unsigned short k, Bdd *bdd, long * head){
	long s = bdd->nodesCount, S = bdd->nodesCount;
	unsigned short j = k, n = bdd->orderTableSize, minpos = k;

	if (k > n / 2){
		for (; j < n - 1; j++){
			S -= swappingVariables(j, bdd, head);
			if ((double)S / s >= 1.25){ break; }
			if (S < s){ s = S; minpos = j + 1; }
		}
		if (j == n - 1){ j--; }
		for (; j > 0; j--){
			S -= swappingVariables(j, bdd, head);
			if ((double)S / s >= 1.25){ break; }
			if (S < s){ s = S; minpos = j; }
		}
		if (j == 0){ j = 1; }
		if (s == S){
			for (; j < minpos; j++){
				S -= swappingVariables(j, bdd, head);
			}
		}
		else{
			for (; S != s; j++){
				S -= swappingVariables(j, bdd, head);
			}
			//j--;
		}
	}
	else{
		for (; j > 1; j--){
			S -= swappingVariables(j - 1, bdd, head);
			if ((double)S / s >= 1.25){ break; }
			if (S < s){ s = S; minpos = j; }
		}
		//if (j == -1){ j = 0; }
		for (; j < n - 1; j++){
			S -= swappingVariables(j, bdd, head);
			if ((double)S / s >= 1.25){ break; }
			if (S < s){ s = S; minpos = j + 1; }
		}
		if (j == n - 1){ j--; }
		if (s == S){
			for (; j >= minpos; j--){
				S -= swappingVariables(j, bdd, head);
			}
		}
		else{
			for (; (S != s) && (j > 0); j--){
				S -= swappingVariables(j, bdd, head);
			}
		}
		j++;
	}
	return j;
}

//reorders the BDD according to the specified order
//returns quantity of the reordered variables +1
//not to use between calls initializeTraversal() and deinitializeTraversal()
long makeBddOrder(Bdd *bdd, unsigned short * order, long orderSize){
	unsigned short i, j, pos = 1;

	long * head = initializeTraversal(bdd);

	for (i = 1; i < orderSize; i++){
		for (j = 1; j < bdd->orderTableSize; j++){
			if (bdd->orderTable[j] == order[i]){
				for (--j; j >= pos; j--){ swappingVariables(j, bdd, head); }
				pos++; break;
			}
		}
	}

	deinitializeTraversal(bdd, head);
	return pos;
}

long compareBdd(Bdd *L, Bdd *R){
	long result = 0;
	long i, j, k, k2, p, p2;
	unsigned short *reservOrder;
	long *headL, *headR;
	Bdd *reordered;
	long reorderedSize;

	if (L->orderTableSize == R->orderTableSize){
		if (L->orderTableSize == 0){
			if (node_(L->root, L).hi == node_(R->root, R).hi){
				result = 1;
			}
		}
		else{
			for (i = 1; i < L->orderTableSize; i++){
				result = 0;
				for (j = 1; j < L->orderTableSize; j++){
					if (L->orderTable[i] == R->orderTable[j]){
						result = 1; break;
					}
				}
				if (result == 0){ break; }
			}
			if (result){
				if (L->nodesCount>R->nodesCount){
					reordered = L; reorderedSize = L->nodesCount;
				}
				else{ reordered = R; reorderedSize = R->nodesCount; }
				reservOrder = (unsigned short*)malloc(sizeof(short)* L->orderTableSize);;
				for (i = 1; i < L->orderTableSize; i++){
					reservOrder[i] = reordered->orderTable[i];
				}
				if (reordered == L){ makeBddOrder(L, R->orderTable, L->orderTableSize); }
				else{ makeBddOrder(R, L->orderTable, L->orderTableSize); }

				if (L->nodesCount != R->nodesCount){ result = 0; }
				else{

					headL = initializeTraversal(L);
					headR = initializeTraversal(R);

					node_(0, L).aux = node_(0, L).lo; node_(1, L).aux = node_(1, L).lo;
					node_(0, R).aux = node_(0, R).lo; node_(1, R).aux = node_(1, R).lo;

					for (i = L->orderTableSize - 1; i > 0; i--){
						p = ~headL[i];
						while (p){
							p2 = ~node_(p, L).aux;
							k = ~headR[i]; k2 = 0; result = 0;
							while (k){
								if ((node_(node_(p, L).lo, L).aux == node_(node_(p, R).lo, R).aux) && (node_(node_(p, L).hi, L).aux == node_(node_(k, R).hi, R).aux)){
									result = 1;
									if (k2 == 0){ headR[i] = node_(k, R).aux; }
									else{ node_(k2, R).aux = node_(k, R).aux; }
									node_(k, R).aux = p;
									break;
								}
								k2 = k; k = ~node_(k, R).aux;
							}
							if (result == 0){
								for (j = i; j > 0; j--){
									k = ~headR[j];
									while (k){
										k2 = ~node_(k, R).aux; node_(k, R).aux = 0; k = k2;
									}
								}
								headL[i] = ~p;
								for (j = i; j > 0; j--){
									k = ~headL[j];
									while (k){
										k2 = ~node_(k, L).aux; node_(k, L).aux = 0; k = k2;
									}
								}
								i = 0; break;
							}
							node_(p, L).aux = p;
							p = p2;
						}
					}

					free(headL); free(headR);
				}
				if (reordered->nodesCount>reorderedSize){ makeBddOrder(reordered, reservOrder, reordered->orderTableSize); }
				free(reservOrder);
			}
		}
	}

	return result;
}

//puts BDD nodes sequentially
//not to use between calls initializeTraversal() and deinitializeTraversal()
void compactBdd(Bdd *bdd){
	BddNode tempNode;
	long s, p, j, currPos, tempAux;
	unsigned long i;

	long * head = (long*)malloc(bdd->orderTableSize*sizeof(long));
	long * lcount = (long*)malloc(bdd->orderTableSize*sizeof(long));

	for (p = 0; p < bdd->orderTableSize; p++) { head[p] = -1; lcount[p] = 0; }
	s = bdd->root; node_(bdd->root, bdd).aux = ~0; node_(bdd->terminals[0], bdd).aux = ~0; node_(bdd->terminals[1], bdd).aux = ~0;
	while (s != 0){
		p = s;
		s = ~node_(p, bdd).aux;
		i = ~head[node_(p, bdd).v]; j = 0;
		if (head[node_(p, bdd).v] != -1){
			while (((long)i < p) && (i != 0)){ j = i; i = ~node_(i, bdd).aux; }
		}
		if (j == 0){
			node_(p, bdd).aux = head[node_(p, bdd).v]; head[node_(p, bdd).v] = ~p;
		}
		else{
			node_(p, bdd).aux = node_(j, bdd).aux; node_(j, bdd).aux = ~p;
		}
		lcount[node_(p, bdd).v]++; node_(p, bdd).ref = 0;
		if (node_(node_(p, bdd).lo, bdd).aux >= 0){ node_(node_(p, bdd).lo, bdd).aux = ~s; s = node_(p, bdd).lo; }
		if (node_(node_(p, bdd).hi, bdd).aux >= 0){ node_(node_(p, bdd).hi, bdd).aux = ~s; s = node_(p, bdd).hi; }
	}

	for (p = bdd->orderTableSize - 2; p > 0; p--){
		lcount[p] += lcount[p + 1];
	}

	for (i = 0; i < 2; i++){
		node_(bdd->terminals[i], bdd).aux = node_(bdd->terminals[i], bdd).lo;
		node_(bdd->terminals[i], bdd).lo = i;
		node_(bdd->terminals[i], bdd).hi = i;
	}

	for (i = 0; i<(bdd->nodesCount) >> LOG_MEMPAGE_SIZE; i++){
		if (bdd->container.page[i] == NULL){ createPage(bdd); }
		bdd->container.page[i]->availableNode = -1;
		bdd->container.page[i]->availableNodeCount = 0;
	}
	s = (bdd->nodesCount - 1) >> LOG_MEMPAGE_SIZE;
	p = s*MEMPAGE_SIZE; bdd->container.page[s]->availableNodeCount = 0;
	for (j = -1, i = (MEMPAGE_SIZE - 1); i>((bdd->nodesCount - 1)&(MEMPAGE_SIZE - 1)); i--){
		if (bdd->container.page[s]->node[i].v == 0){
			bdd->container.page[s]->node[i].hi = j; j = p + i; bdd->container.page[s]->availableNodeCount++;
		}
	}
	bdd->container.page[s]->availableNode = j;

	currPos = 2;
	for (i = bdd->orderTableSize - 1; i > 0; i--){
		j = ~head[i];
		while (j){
			if (j != currPos){
				if (node_(currPos, bdd).aux >= 0){
					node_(currPos, bdd).v = node_(j, bdd).v;
				}
				else{
					tempNode.v = node_(currPos, bdd).v;	tempNode.lo = node_(currPos, bdd).lo;
					tempNode.hi = node_(currPos, bdd).hi; tempNode.aux = node_(currPos, bdd).aux;

					tempAux = 2 + lcount[node_(currPos, bdd).v + 1];
					s = ~head[node_(currPos, bdd).v]; p = 0;
					while (s != currPos){ p = s; s = ~node_(s, bdd).aux; tempAux++; }
					node_(currPos, bdd).aux = tempAux;
					if (node_(tempAux, bdd).aux != 0){
						tempAux = allocateNode(bdd);
						node_(tempAux, bdd).ref = 1;
					}
					if (p == 0){ head[node_(currPos, bdd).v] = ~tempAux; }
					else { node_(p, bdd).aux = ~tempAux; }

					node_(currPos, bdd).v = node_(j, bdd).v;
					node_(tempAux, bdd).v = tempNode.v; node_(tempAux, bdd).lo = tempNode.lo;
					node_(tempAux, bdd).hi = tempNode.hi; node_(tempAux, bdd).aux = tempNode.aux;
				}
			}
			s = j;
			node_(currPos, bdd).lo = node_(node_(j, bdd).lo, bdd).aux;
			node_(currPos, bdd).hi = node_(node_(j, bdd).hi, bdd).aux;
			j = ~node_(j, bdd).aux;
			head[i] = ~j;
			if (node_(s, bdd).ref){
				node_(s, bdd).ref = 0; freeNode(s, bdd);
			}
			else{
				node_(s, bdd).aux = currPos;
			}
			currPos++;
		}
	}

	bdd->root = bdd->nodesCount - 1;
	free(head);

	s = (bdd->nodesCount - 1) >> LOG_MEMPAGE_SIZE;
	p = s*MEMPAGE_SIZE;
	for (i = s + 1; i < (unsigned long)bdd->container.pageCount; i++){
		free(bdd->container.page[i]);
		bdd->container.page[i] = NULL;
	}
	j = bdd->container.page[s]->availableNode;
	for (i = (MEMPAGE_SIZE - 1); i>((bdd->nodesCount - 1)&(MEMPAGE_SIZE - 1)); i--){
		if (bdd->container.page[s]->node[i].v != 0){
			bdd->container.page[s]->node[i].hi = j; j = p + i; bdd->container.page[s]->availableNodeCount++;
			bdd->container.page[s]->node[i].v = 0;
		}
	}
	bdd->container.page[s]->availableNode = j;
}

//binaryMode : 0 - text, 1 - binary. 
//the program doesn't read the text (for a BDD reuse use only binary data)
//*head can be useful for debugging of some functions
void fprintBdd(FILE * file, Bdd *bdd, long binaryMode){
	long i, j;
	long f, g, p, s;
	unsigned short index;
	long *head;

	if (!binaryMode){
		fprintf(file, "BDD size: %d\n", bdd->nodesCount);
		fprintf(file, "nid - v ? lo : hi\n");

		if ((bdd->root == bdd->terminals[0]) || (bdd->root == bdd->terminals[1])){
			fprintf(file, "%d - 0 ? %d : %d\n", bdd->root, node_(bdd->root, bdd).lo, node_(bdd->root, bdd).hi);
		}
		else{
			head = (long*)malloc(bdd->orderTableSize*sizeof(long));
			for (p = 0; p<bdd->orderTableSize; p++) { head[p] = -1; }
			s = bdd->root; node_(bdd->root, bdd).aux = ~0; node_(bdd->terminals[0], bdd).aux = ~0; node_(bdd->terminals[1], bdd).aux = ~0;
			while (s != 0){
				p = s;
				s = ~node_(p, bdd).aux;
				i = ~head[node_(p, bdd).v]; j = 0;
				if (head[node_(p, bdd).v] != -1){
					while (i > p){ j = i; i = ~node_(i, bdd).aux; }
				}
				if (j == 0){
					node_(p, bdd).aux = head[node_(p, bdd).v]; head[node_(p, bdd).v] = ~p;
				}
				else{
					node_(p, bdd).aux = node_(j, bdd).aux; node_(j, bdd).aux = ~p;
				}
				if (node_(node_(p, bdd).lo, bdd).aux >= 0){ node_(node_(p, bdd).lo, bdd).aux = ~s; s = node_(p, bdd).lo; }
				if (node_(node_(p, bdd).hi, bdd).aux >= 0){ node_(node_(p, bdd).hi, bdd).aux = ~s; s = node_(p, bdd).hi; }
			}
			for (i = 1; i < bdd->orderTableSize; i++){
				j = ~head[i];
				while (j){
					f = node_(j, bdd).lo; if (f < 2){ f = node_(f, bdd).lo; }
					g = node_(j, bdd).hi; if (g < 2){ g = node_(g, bdd).lo; }
					index = bdd->orderTable[node_(j, bdd).v];
					fprintf(file, "%d - %c%d ? %d : %d\n", j, varletter_(index), varindex_(index), f, g);
					j = ~node_(j, bdd).aux;
				}
			}
			deinitializeTraversal(bdd, head);
		}
	}
	else{
		compactBdd(bdd);
		fwrite(&(bdd->orderTableSize), sizeof(short), 1, file);
		for (i = 1; i < bdd->orderTableSize; i++){
			fwrite(&(bdd->orderTable[i]), sizeof(short), 1, file);
		}
		fwrite(&(bdd->nodesCount), sizeof(long), 1, file);
		if (bdd->root<2){ fwrite(&(node_(bdd->root, bdd).lo), sizeof(long), 1, file); }
		else{
			fwrite(&(bdd->root), sizeof(long), 1, file);
			for (i = 2; i < bdd->nodesCount; i++){
				f = node_(i, bdd).lo; if (f < 2){ f = node_(f, bdd).lo; }
				g = node_(i, bdd).hi; if (g < 2){ g = node_(g, bdd).lo; }
				fwrite(&(node_(i, bdd).v), sizeof(short), 1, file);
				fwrite(&f, sizeof(long), 1, file);
				fwrite(&g, sizeof(long), 1, file);
			}
		}
	}
}

Bdd* loadBdd(FILE* file){
	long i;

	Bdd* bdd = (Bdd*)calloc(1, sizeof(Bdd));
	allocateNode(bdd); node_(0, bdd).v = -1;
	node_(0, bdd).lo = 0; node_(0, bdd).hi = 0;
	allocateNode(bdd); node_(1, bdd).v = -1;
	node_(1, bdd).lo = 1; node_(1, bdd).hi = 1;
	bdd->terminals[0] = 0; bdd->terminals[1] = 1;

	fread(&(bdd->orderTableSize), sizeof(short), 1, file);
	bdd->orderTable = (unsigned short*)malloc(sizeof(short)*bdd->orderTableSize);
	for (i = 1; i < bdd->orderTableSize; i++){
		fread(&(bdd->orderTable[i]), sizeof(short), 1, file);
	}
	fread(&(bdd->nodesCount), sizeof(long), 1, file);
	fread(&(bdd->root), sizeof(long), 1, file);
	for (i = 2; i < bdd->nodesCount; i++){
		allocateNode(bdd);
		fread(&(node_(i, bdd).v), sizeof(short), 1, file);
		fread(&(node_(i, bdd).lo), sizeof(long), 1, file);
		fread(&(node_(i, bdd).hi), sizeof(long), 1, file);
	}
	return bdd;
}

//search of an optimum order of variables
//returns quantity of the disposed nodes
//not to use between calls initializeTraversal() and deinitializeTraversal()
long optimizeBdd(Bdd *bdd){
	long *head;	unsigned short * reorderTable;

	unsigned short i, j, k, oldPos, newPos;
	long startSize = bdd->nodesCount;

	head = initializeTraversal(bdd);
	reorderTable = (unsigned short*)malloc(sizeof(short)*bdd->orderTableSize);
	for (i = 1; i < bdd->orderTableSize; i++){ reorderTable[i] = i; }
	srand(time(NULL));

	for (k = 0; k < 1; k++){
		for (i = bdd->orderTableSize - 1; i > 0; i--){
			oldPos = (rand() % i) + 1;
			j = reorderTable[i]; reorderTable[i] = reorderTable[oldPos];
			reorderTable[oldPos] = j;
		}
		for (i = 1; i < bdd->orderTableSize; i++){
			//for (i = bdd.orderTableSize-1; i>=0; i--){
			for (j = 1; j < bdd->orderTableSize; j++){
				if (reorderTable[j] == i){ oldPos = j; break; }
			}

			newPos = siftingVariable(oldPos, bdd, head);
			if (newPos < oldPos){
				for (j = oldPos; j>newPos; j--){ reorderTable[j] = reorderTable[j - 1]; }
				reorderTable[newPos] = i;
			}
			else if (newPos>oldPos){
				for (j = oldPos; j < newPos; j++){ reorderTable[j] = reorderTable[j + 1]; }
				reorderTable[newPos] = i;
			}
		}
	}
	free(reorderTable);
	deinitializeTraversal(bdd, head);

	compactBdd(bdd);

	return startSize - bdd->nodesCount;
}

long copyBdd(Bdd *target, Bdd *source){
	long currPos, j, k;
	long * head = initializeTraversal(source);
	unsigned short i;

	for (i = 0; i < target->container.pageCount; i++){
		if (target->container.page[i] != NULL){
			free(target->container.page[i]);
		}
	}
	free(target->container.page); target->container.pageCount = 0;
	free(target->orderTable);
	target->orderTable = (unsigned short*)malloc(sizeof(short)*source->orderTableSize);
	for (i = 1; i < source->orderTableSize; i++){
		target->orderTable[i] = source->orderTable[i];
	}
	target->orderTableSize = source->orderTableSize;

	for (i = 0; i < 2; i++){
		currPos = allocateNode(target);
		node_(currPos, target).v = -1;
		node_(currPos, target).lo = i;
		node_(currPos, target).hi = i;
		node_(source->terminals[i], source).aux = currPos;
	}

	for (i = source->orderTableSize - 1; i > 0; i--){
		j = ~head[i];
		while (j){
			currPos = allocateNode(target);
			node_(currPos, target).v = i;
			node_(currPos, target).lo = node_(node_(j, source).lo, source).aux;
			node_(currPos, target).hi = node_(node_(j, source).hi, source).aux;;
			k = j;
			j = ~node_(j, source).aux;
			node_(k, source).aux = currPos;
		}
	}

	target->nodesCount = source->nodesCount;
	target->root = currPos;
	free(head);//deinitializeTraversal is impossible, head table is corrupted
	return currPos + 1;
}

long findLevel(long f, long g, char op, Bdd *L, Bdd *R){
	if (f <= 1){ f = node_(f, L).lo; }
	if (g <= 1){ g = node_(g, R).lo; }
	long t;
	if ((f <= 1) && (g <= 1)){ return -((op >> (3 - 2 * node_(f, L).lo - node_(g, R).lo)) & 1); }
	if ((f <= 1) && (g > 1)){
		t = (node_(f, L).lo ? op & 3 : op >> 2);
		if (t == 0){ return 0; }
		if (t == 3){ return -1; }
	}
	if ((f > 1) && (g <= 1)){
		t = (node_(g, R).lo ? op : op >> 1) & 5;
		if (t == 0){ return 0; }
		if (t == 5){ return -1; }
	}
	return (node_(f, L).v>node_(g, R).v ? node_(g, R).v : node_(f, L).v);
}

long makeTemplate(long f, long g, long *tbot, long hbase, long logTableSize, TemplateContainer *container, Bdd* L, Bdd *R){
	long t, h;
	h = hbase - getHash(f, g, logTableSize);
	t = template_(h, *container).h;
	while ((t != -1) && ((template_(t, *container).left != f) || (template_(t, *container).right != g))){
		t = template_(t, *container).l;
	}
	if (t == -1){
		t = (*tbot) + 1; (*tbot) = t;
		template_(t, *container).left = f;
		template_(t, *container).right = g;
		template_(t, *container).l = template_(h, *container).h;
		template_(h, *container).h = t;
	}
	return t;
}

//creates the BDD L op R (L+R, L*R or other)
//op - operation truth table.
//				L	0011
//				R	0101
//Sample: AND - 1  (0001); OR - 7 (0111); XOR - 6 (0110)
Bdd* synthesizeBdd(Bdd *L, Bdd *R, char op, long saveL, long saveR){
	Bdd *bdd = (Bdd*)calloc(1, sizeof(Bdd));

	unsigned short * orderBackup = NULL, *pseudoOrder = NULL;

	long * head;
	long i, j, k = 0;
	Bdd *smaller, *bigger; long saveSmaller;

	//initiation
	allocateNode(bdd); node_(0, bdd).v = -1;
	node_(0, bdd).lo = 0; node_(0, bdd).hi = 0;
	allocateNode(bdd); node_(1, bdd).v = -1;
	node_(1, bdd).lo = 1; node_(1, bdd).hi = 1;
	bdd->terminals[0] = 0; bdd->terminals[1] = 1;
	bdd->nodesCount = 2;

	if ((L->root < 2) || (R->root < 2)){
		if ((L->root < 2) && (R->root < 2)){
			bdd->root = ((op >> (3 - 2 * node_(L->root, L).lo - node_(R->root, R).lo)) & 1);
		}
		else if (L->root<2){
			op = (node_(L->root, L).lo ? op & 3 : op >> 2);
			if (op == 0){ bdd->root = 0; }
			else if (op == 3){ bdd->root = 1; }
			else {
				if (saveR){ copyBdd(bdd, R); }
				else{
					free(bdd->orderTable);
					free(bdd->container.page[0]);
					free(bdd->container.page);
					bdd->container.page = R->container.page;
					bdd->container.pageCount = R->container.pageCount;
					bdd->orderTable = R->orderTable;
					bdd->orderTableSize = R->orderTableSize;
					bdd->nodesCount = R->nodesCount;
					bdd->root = R->root;
					bdd->terminals[0] = R->terminals[0];
					bdd->terminals[1] = R->terminals[1];
					R->container.page = NULL; R->container.pageCount = 0;
					R->orderTable = NULL;
				}
				if (op == 2){ inverseBdd(bdd); }
			}
		}
		else{
			op = (node_(R->root, R).lo ? op : op >> 1) & 5;
			if (op == 0){ bdd->root = 0; }
			else if (op == 5){ bdd->root = 1; }
			else {
				if (saveL){ copyBdd(bdd, L); }
				else{
					free(bdd->orderTable);
					free(bdd->container.page[0]);
					free(bdd->container.page);
					bdd->container.page = L->container.page;
					bdd->container.pageCount = L->container.pageCount;
					bdd->orderTable = L->orderTable;
					bdd->orderTableSize = L->orderTableSize;
					bdd->nodesCount = L->nodesCount;
					bdd->root = L->root;
					bdd->terminals[0] = L->terminals[0];
					bdd->terminals[1] = L->terminals[1];
					L->container.page = NULL; L->container.pageCount = 0;
					L->orderTable = NULL;
				}
				if (op == 4){ inverseBdd(bdd); }
			}
		}
		if (!saveL){ freeBdd(L); }
		if (!saveR){ freeBdd(R); }
	}
	else{
		if (L->orderTableSize >= R->orderTableSize){
			smaller = R; bigger = L; saveSmaller = saveR;
		}
		else{
			smaller = L; bigger = R; saveSmaller = saveL;
		}

		if (saveSmaller){
			orderBackup = (unsigned short*)malloc(sizeof(short)* smaller->orderTableSize);
			for (i = 1; i < smaller->orderTableSize; i++){ orderBackup[i] = smaller->orderTable[i]; }
		}
		k = -1;
		k = makeBddOrder(smaller, bigger->orderTable, bigger->orderTableSize);
		bdd->orderTableSize = L->orderTableSize + (R->orderTableSize - 1) - (k - 1);
		bdd->orderTable = (unsigned short*)malloc(sizeof(short)*bdd->orderTableSize);
		for (i = 1; i < bigger->orderTableSize; i++){
			bdd->orderTable[i] = bigger->orderTable[i];
		}
		pseudoOrder = (unsigned short*)malloc(sizeof(short)* smaller->orderTableSize);
		for (j = 1, i = 1; i < k; i++){
			for (; j < bigger->orderTableSize; j++){
				if (smaller->orderTable[i] == bigger->orderTable[j]){
					pseudoOrder[i] = j; break;
				}
			}
		}
		for (; i < smaller->orderTableSize; i++){
			pseudoOrder[i] = bigger->orderTableSize + i - k;
			bdd->orderTable[bigger->orderTableSize + i - k] = smaller->orderTable[i];
		}


		head = initializeTraversal(smaller);
		for (i = 1; i < smaller->orderTableSize; i++){
			j = ~head[i];
			while (j){
				node_(j, smaller).v = pseudoOrder[i];
				j = ~node_(j, smaller).aux;
			}
		}
		free(pseudoOrder);

		//algorithm S (chapter 7.1.4)

		long *lstart, *lcount, *llist, *hlist;
		long f, g, l, ll, lh, k, t, s, r, tbot, ntop, hbase, b;
		long lnode, rnode;
		unsigned short vf, vg;

		lstart = (long*)malloc(sizeof(long)*bdd->orderTableSize);
		lcount = (long*)malloc(sizeof(long)*bdd->orderTableSize);
		llist = (long*)malloc(sizeof(long)*bdd->orderTableSize);
		hlist = (long*)malloc(sizeof(long)*bdd->orderTableSize);

		TemplateContainer container = { 0 }; allocateTemplate(&container, 1);
		template_(0, container).left = ~0;
		template_(0, container).right = 0;
		template_(1, container).left = ~1;
		template_(1, container).right = 1;

		//S1
		f = L->root; g = R->root; l = findLevel(f, g, op, L, R);
		lstart[l - 1] = 1;
		for (k = l; k < bdd->orderTableSize; k++){
			llist[k] = -1; hlist[k] = -1; lcount[k] = 0;
		}
		tbot = 2; allocateTemplate(&container, 2);
		template_(tbot, container).left = f; template_(tbot, container).right = g;

		do{
			i = 0;
			do{
				//S2
				lstart[l] = tbot; t = lstart[l - 1];
				while (t < tbot){
					t++; allocateTemplate(&container, t);
					f = template_(t, container).left; g = template_(t, container).right;
					vf = node_(f, L).v; vg = node_(g, R).v;

					lnode = ((vf <= vg) && (f > 1) ? (node_(f, L).lo > 1 ? node_(f, L).lo : node_(node_(f, L).lo, L).lo) : f);
					rnode = ((vf >= vg) && (g > 1) ? (node_(g, R).lo > 1 ? node_(g, R).lo : node_(node_(g, R).lo, R).lo) : g);
					ll = findLevel(lnode, rnode, op, L, R);
					lnode = ((vf <= vg) && (f > 1) ? (node_(f, L).hi > 1 ? node_(f, L).hi : node_(node_(f, L).hi, L).hi) : f);
					rnode = ((vf >= vg) && (g > 1) ? (node_(g, R).hi > 1 ? node_(g, R).hi : node_(node_(g, R).hi, R).hi) : g);
					lh = findLevel(lnode, rnode, op, L, R);

					if (ll <= 0){ template_(t, container).l = -ll; }
					else{ template_(t, container).l = llist[ll]; llist[ll] = t; lcount[ll]++; }
					if (lh <= 0){ template_(t, container).h = -lh; }
					else{ template_(t, container).h = hlist[lh]; hlist[lh] = t; lcount[lh]++; }
				}
				//S3
				if (l == bdd->orderTableSize - 1){
					i = 1; break;
				}
				l++;
			} while (lcount[l] == 0);
			if (i && (l == bdd->orderTableSize - 1)){ break; }
			//S4
			for (hbase = 2, b = 1; hbase < lcount[l]; hbase <<= 1, b++);
			hbase += tbot; allocateTemplate(&container, hbase);
			for (k = tbot + 1; k <= hbase; k++){ template_(k, container).h = -1; }
			//S5
			t = llist[l];
			while (t != -1){
				s = template_(t, container).l;
				f = template_(t, container).left;
				g = template_(t, container).right;
				vf = node_(f, L).v;
				vg = node_(g, R).v;
				lnode = ((vf <= vg) && (f > 1) ? (node_(f, L).lo > 1 ? node_(f, L).lo : node_(node_(f, L).lo, L).lo) : f);
				rnode = ((vf >= vg) && (g > 1) ? (node_(g, R).lo > 1 ? node_(g, R).lo : node_(node_(g, R).lo, R).lo) : g);
				template_(t, container).l = makeTemplate(lnode, rnode, &tbot, hbase, b, &container, L, R);
				t = s;
			}
			t = hlist[l];
			while (t != -1){
				s = template_(t, container).h;
				f = template_(t, container).left;
				g = template_(t, container).right;
				vf = node_(f, L).v;
				vg = node_(g, R).v;
				lnode = ((vf <= vg) && (f > 1) ? (node_(f, L).hi > 1 ? node_(f, L).hi : node_(node_(f, L).hi, L).hi) : f);
				rnode = ((vf >= vg) && (g > 1) ? (node_(g, R).hi > 1 ? node_(g, R).hi : node_(node_(g, R).hi, R).hi) : g);
				template_(t, container).h = makeTemplate(lnode, rnode, &tbot, hbase, b, &container, L, R);
				t = s;
			}
		} while (1);
		//S6
		ntop = 2;
		g = 0;
		//{deleting operands
		if (saveSmaller){
			for (i = 1; i < smaller->orderTableSize; i++){
				j = ~head[i];
				while (j){
					node_(j, smaller).v = (unsigned short)i;
					j = ~node_(j, smaller).aux;
				}
			}
			deinitializeTraversal(smaller, head);
			makeBddOrder(smaller, orderBackup, smaller->orderTableSize);
			free(orderBackup);
		}
		else{
			freeBdd(smaller); free(head);
		}
		if (bigger == L){ if (!saveL){ freeBdd(L); } }
		if (bigger == R){ if (!saveR){ freeBdd(R); } }
		//}
		do{
			//S7
			t = lstart[l - 1];
			while (t < lstart[l]){
				t++;
				template_(t, container).l = template_(template_(t, container).l, container).right;
				template_(t, container).h = template_(template_(t, container).h, container).right;
				if (template_(t, container).l == template_(t, container).h){ template_(t, container).right = template_(t, container).l; }
				else{
					template_(t, container).right = -1;
					template_(t, container).left = template_(template_(t, container).l, container).left;
					template_(template_(t, container).l, container).left = t;
				}
			}
			//S8
			while (1){
				if (t == lstart[l - 1]){ t = lstart[l] + 1; break; }
				if (template_(t, container).left < 0){
					template_(template_(t, container).l, container).left = template_(t, container).left;
				}
				t--;
			}
			while (1){
				//S9
				i = 0;
				do{
					t--;
					if (t == lstart[l - 1]) { i = 1; break; }
				} while (template_(t, container).right >= 0);
				if ((t == lstart[l - 1]) && i) { break; }
				//S10
				s = t;
				while (s > 0){
					r = template_(s, container).h;
					template_(s, container).right = template_(r, container).left;
					if (template_(r, container).left < 0){ template_(r, container).left = s; }
					s = template_(s, container).left;
				}
				s = t;
				//S11
				while (s >= 0){
					if (template_(s, container).right >= 0){ s = template_(s, container).left; }
					else{
						r = template_(s, container).left;
						template_(template_(s, container).h, container).left = template_(s, container).right;
						template_(s, container).right = s;
						g = ntop; allocateNode(bdd); ntop++;
						bdd->nodesCount++;
						template_(s, container).left = ~g;
						node_(g, bdd).lo = ~template_(template_(s, container).l, container).left;
						node_(g, bdd).hi = ~template_(template_(s, container).h, container).left;
						node_(g, bdd).v = (unsigned short)l;
						s = r;
					}
				}
			}
			//S12
			l--;
		} while (lstart[l] > 1);
		if (template_(2, container).right == 0){ ntop--; }
		freeTemplates(&container);

		bdd->root = ntop - 1;
		if (ntop < 2){ free(bdd->orderTable); bdd->orderTable = 0; bdd->orderTableSize = 0; }
		checkBddOrderTable(bdd);
	}
	return bdd;
}

#ifdef __windows__
DWORD WINAPI createBddThread(LPVOID lpParam){
	long index = ((BTreeSubContainer*)lpParam)->index;
	BTreeContainer *tree = ((BTreeSubContainer*)lpParam)->tree;
	free(lpParam);
	Bdd *bdd;

	long i, j = index;
	while (treeNode_(j).right != -1){
		j = treeNode_(j).right;
	}

	bdd = (Bdd*)calloc(1, sizeof(Bdd));
	allocateNode(bdd); node_(0, bdd).v = -1;
	node_(0, bdd).lo = 0; node_(0, bdd).hi = 0;
	allocateNode(bdd); node_(1, bdd).v = -1;
	node_(1, bdd).lo = 1; node_(1, bdd).hi = 1;
	bdd->terminals[0] = 0; bdd->terminals[1] = 1;
	if (treeNode_(j).op > '1'){
		bdd->nodesCount = 3; bdd->orderTableSize = 2;
		bdd->orderTable = (unsigned short*)malloc(2 * sizeof(short));
		bdd->orderTable[1] = varhash_(treeNode_(j).op, treeNode_(j).index);
		allocateNode(bdd); node_(2, bdd).v = 1;
		node_(2, bdd).lo = 0; node_(2, bdd).hi = 1;
		bdd->root = 2;
	}
	else{
		bdd->nodesCount = 2; bdd->orderTableSize = 0; bdd->orderTable = NULL;
		bdd->root = treeNode_(j).op - '0';
	}

	treeNode_(j).bdd = bdd; treeNode_(j).left = bdd->nodesCount;

	while (j != index){
		j = treeNode_(j).parent;
		if (treeNode_(j).left != -1){
			WaitForSingleObject(treeNode_(treeNode_(j).left).thread, INFINITE);
			CloseHandle(treeNode_(treeNode_(j).left).thread);
			treeNode_(j).bdd = synthesizeBdd((Bdd*)treeNode_(treeNode_(j).left).bdd, (Bdd*)treeNode_(treeNode_(j).right).bdd, treeNode_(j).op, 0, 0);
			i = treeNode_(j).left;
			treeNode_(j).left = (treeNode_(treeNode_(j).left).left > treeNode_(treeNode_(j).right).left ? treeNode_(treeNode_(j).left).left : treeNode_(treeNode_(j).right).left);
			if ((treeNode_(j).left * 2) < treeNode_(j).bdd->nodesCount){ optimizeBdd((Bdd*)treeNode_(j).bdd); treeNode_(j).left = treeNode_(j).bdd->nodesCount; }
			freeBTreeNode(tree, i);
		}
		else{
			treeNode_(j).bdd = treeNode_(treeNode_(j).right).bdd;
			inverseBdd((Bdd*)treeNode_(j).bdd);
		}
		freeBTreeNode(tree, treeNode_(j).right);
	}
	return 0;
}

Bdd* createBdd(char *boolFunc){
	Bdd *bdd = NULL;
	char operators[] = { ' ', '&', '>', '!', '<', '@', '^', '|', '_', 'A', 'N', '#', 'Y', 'S', 'D', 'E' };

	long i, j = 0, k = 0;
	BTreeContainer *tree = (BTreeContainer *)calloc(1, sizeof(BTreeContainer));
	InitializeCriticalSection(&(tree->critSec));
	tree->node = (BTreeNode*)malloc(1024 * sizeof(BTreeNode));
	tree->availableNode = -1;
	for (i = 1023; i > 0; i--){
		treeNode_(i).parent = tree->availableNode;
		tree->availableNode = i;
		treeNode_(i).bdd = NULL;
		treeNode_(i).thread = NULL;
	}
	treeNode_(0).index = 0;
	BTreeSubContainer *subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
	subTree->tree = tree;

	long tempPos, currStringPos = 0, currTreePos = 1, oldTreePos = 0, offset = 0;
	char currS;

	DWORD threadId;

	//parse of request
	currS = boolFunc[currStringPos];
	while (currS>'\n'){
		currTreePos = allocateBTreeNode(tree);
		if (j == 0){
			treeNode_(oldTreePos).right = currTreePos;
			treeNode_(currTreePos).parent = oldTreePos;
			treeNode_(currTreePos).left = -1;
			if (currS == '~'){
				treeNode_(currTreePos).op = 0;
				treeNode_(currTreePos).index = offset + 3;

			}
			else if ((currS == '(') || (currS == '[') || (currS == '{')){
				offset += 4; currStringPos++;
				currS = boolFunc[currStringPos];
				continue;
			}
			else{
				if (((currS >= '0') && (currS <= '1')) || ((currS >= 'A') && (currS <= 'Z')) || ((currS >= 'a') && (currS <= 'z'))){
					treeNode_(currTreePos).right = -1;
					treeNode_(currTreePos).op = currS;
					treeNode_(currTreePos).index = 0;
					currStringPos++;
					while ((boolFunc[currStringPos] >= '0') && (boolFunc[currStringPos] <= '9')){
						treeNode_(currTreePos).index = treeNode_(currTreePos).index * 10 + (boolFunc[currStringPos] - '0');
						currStringPos++;
					}
					currStringPos--;
					j = 1;
				}
				else{ break; }
			}
		}
		else if (j == 1){
			if ((currS == ')') || (currS == ']') || (currS == '}')){
				offset -= 4; currStringPos++;
				currS = boolFunc[currStringPos];
				if (offset < 0){ j = 0; break; }
				continue;
			}
			else{
				j = 0;
				for (i = 1; i < 16; i++){ if (operators[i] == currS){ break; } }
				if (i < 16){
					treeNode_(currTreePos).op = (char)i;
					if (i == 1){ treeNode_(currTreePos).index = offset + 3; }
					else if ((i == 6) || (i == 7)){ treeNode_(currTreePos).index = offset + 2; }
					else { treeNode_(currTreePos).index = offset + 1; }
					oldTreePos = treeNode_(oldTreePos).parent;
					while (treeNode_(oldTreePos).index >= treeNode_(currTreePos).index){
						oldTreePos = treeNode_(oldTreePos).parent;
					}
					tempPos = treeNode_(oldTreePos).right;
					treeNode_(oldTreePos).right = currTreePos;
					treeNode_(tempPos).parent = currTreePos;
					treeNode_(currTreePos).parent = oldTreePos;
					treeNode_(currTreePos).left = tempPos;
					subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
					subTree->tree = tree;
					subTree->index = tempPos;
					treeNode_(tempPos).thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)createBddThread, subTree, 0, &threadId);
					if (treeNode_(tempPos).thread == NULL){ ExitProcess(3); }
				}
				else{ break; }
			}
		}
		currStringPos++; currS = boolFunc[currStringPos];
		oldTreePos = currTreePos;
	}

	if (j == 0){ //incorrect request
		for (i = 1; i < 1024; i++){
			if (treeNode_(i).thread != NULL){
				WaitForSingleObject(treeNode_(i).thread, INFINITE);
			}
		}
		for (i = 1; i < 1024; i++){
			if (treeNode_(i).thread != NULL){
				CloseHandle(treeNode_(i).thread);
			}
			if (treeNode_(i).bdd != NULL){
				freeBdd((Bdd*)treeNode_(i).bdd);
				treeNode_(i).bdd = NULL;
			}
		}
	}
	else{
		i = treeNode_(0).right;
		subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
		subTree->tree = tree;
		subTree->index = i;
		treeNode_(i).thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)createBddThread, subTree, 0, &threadId);
		WaitForSingleObject(treeNode_(i).thread, INFINITE);
		CloseHandle(treeNode_(i).thread);
	}

	bdd = (Bdd*)treeNode_(treeNode_(0).right).bdd;
	DeleteCriticalSection(&(tree->critSec));
	free(tree->node); free(tree);
	return bdd;
}
#endif

#ifdef __linux__
void* createBddThread(void *lpParam){
	long index = ((BTreeSubContainer*)lpParam)->index;
	BTreeContainer *tree = ((BTreeSubContainer*)lpParam)->tree;
	free(lpParam);
	Bdd *bdd;

	long i, j = index;
	while (treeNode_(j).right != -1){
		j = treeNode_(j).right;
	}

	bdd = (Bdd*)calloc(1, sizeof(Bdd));
	allocateNode(bdd); node_(0, bdd).v = -1;
	node_(0, bdd).lo = 0; node_(0, bdd).hi = 0;
	allocateNode(bdd); node_(1, bdd).v = -1;
	node_(1, bdd).lo = 1; node_(1, bdd).hi = 1;
	bdd->terminals[0] = 0; bdd->terminals[1] = 1;
	if (treeNode_(j).op > '1'){
		bdd->nodesCount = 3; bdd->orderTableSize = 2;
		bdd->orderTable = (unsigned short*)malloc(2 * sizeof(short));
		bdd->orderTable[1] = varhash_(treeNode_(j).op, treeNode_(j).index);
		allocateNode(bdd); node_(2, bdd).v = 1;
		node_(2, bdd).lo = 0; node_(2, bdd).hi = 1;
		bdd->root = 2;
	}
	else{
		bdd->nodesCount = 2; bdd->orderTableSize = 0; bdd->orderTable = NULL;
		bdd->root = treeNode_(j).op - '0';
	}

	treeNode_(j).bdd = bdd; treeNode_(j).left = bdd->nodesCount;

	while (j != index){
		j = treeNode_(j).parent;
		if (treeNode_(j).left != -1){
			pthread_join(treeNode_(treeNode_(j).left).thread, NULL);
			//CloseHandle(treeNode_(treeNode_(j).left).thread);
			treeNode_(j).bdd = synthesizeBdd((Bdd*)treeNode_(treeNode_(j).left).bdd, (Bdd*)treeNode_(treeNode_(j).right).bdd, treeNode_(j).op, 0, 0);
			i = treeNode_(j).left;
			treeNode_(j).left = (treeNode_(treeNode_(j).left).left > treeNode_(treeNode_(j).right).left ? treeNode_(treeNode_(j).left).left : treeNode_(treeNode_(j).right).left);
			if ((treeNode_(j).left * 2) < treeNode_(j).bdd->nodesCount){ optimizeBdd((Bdd*)treeNode_(j).bdd); treeNode_(j).left = treeNode_(j).bdd->nodesCount; }
			freeBTreeNode(tree, i);
		}
		else{
			treeNode_(j).bdd = treeNode_(treeNode_(j).right).bdd;
			inverseBdd((Bdd*)treeNode_(j).bdd);
		}
		freeBTreeNode(tree, treeNode_(j).right);
	}
	return NULL;
}

Bdd* createBdd(char *boolFunc){
	Bdd *bdd = NULL;
	char operators[] = { ' ', '&', '>', '!', '<', '@', '^', '|', '_', 'A', 'N', '#', 'Y', 'S', 'D', 'E' };

	long i, j = 0, k = 0;
	BTreeContainer *tree = (BTreeContainer *)calloc(1, sizeof(BTreeContainer));
	pthread_mutex_init(&(tree->mutex));
	tree->node = (BTreeNode*)malloc(1024 * sizeof(BTreeNode));
	tree->availableNode = -1;
	for (i = 1023; i > 0; i--){
		treeNode_(i).parent = tree->availableNode;
		tree->availableNode = i;
		treeNode_(i).bdd = NULL;
		treeNode_(i).thread = NULL;
	}
	treeNode_(0).index = 0;
	BTreeSubContainer *subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
	subTree->tree = tree;

	long tempPos, currStringPos = 0, currTreePos = 1, oldTreePos = 0, offset = 0;
	char currS;

	int threadId;

	//parse of request
	currS = boolFunc[currStringPos];
	while (currS>'\n'){
		currTreePos = allocateBTreeNode(tree);
		if (j == 0){
			treeNode_(oldTreePos).right = currTreePos;
			treeNode_(currTreePos).parent = oldTreePos;
			treeNode_(currTreePos).left = -1;
			if (currS == '~'){
				treeNode_(currTreePos).op = 0;
				treeNode_(currTreePos).index = offset + 3;

			}
			else if ((currS == '(') || (currS == '[') || (currS == '{')){
				offset += 4; currStringPos++;
				currS = boolFunc[currStringPos];
				continue;
			}
			else{
				if (((currS >= '0') && (currS <= '1')) || ((currS >= 'A') && (currS <= 'Z')) || ((currS >= 'a') && (currS <= 'z'))){
					treeNode_(currTreePos).right = -1;
					treeNode_(currTreePos).op = currS;
					treeNode_(currTreePos).index = 0;
					currStringPos++;
					while ((boolFunc[currStringPos] >= '0') && (boolFunc[currStringPos] <= '9')){
						treeNode_(currTreePos).index = treeNode_(currTreePos).index * 10 + (boolFunc[currStringPos] - '0');
						currStringPos++;
					}
					currStringPos--;
					j = 1;
				}
				else{ break; }
			}
		}
		else if (j == 1){
			if ((currS == ')') || (currS == ']') || (currS == '}')){
				offset -= 4; currStringPos++;
				currS = boolFunc[currStringPos];
				if (offset < 0){ j = 0; break; }
				continue;
			}
			else{
				j = 0;
				for (i = 1; i < 16; i++){ if (operators[i] == currS){ break; } }
				if (i < 16){
					treeNode_(currTreePos).op = (char)i;
					if (i == 1){ treeNode_(currTreePos).index = offset + 3; }
					else if ((i == 6) || (i == 7)){ treeNode_(currTreePos).index = offset + 2; }
					else { treeNode_(currTreePos).index = offset + 1; }
					oldTreePos = treeNode_(oldTreePos).parent;
					while (treeNode_(oldTreePos).index >= treeNode_(currTreePos).index){
						oldTreePos = treeNode_(oldTreePos).parent;
					}
					tempPos = treeNode_(oldTreePos).right;
					treeNode_(oldTreePos).right = currTreePos;
					treeNode_(tempPos).parent = currTreePos;
					treeNode_(currTreePos).parent = oldTreePos;
					treeNode_(currTreePos).left = tempPos;
					subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
					subTree->tree = tree;
					subTree->index = tempPos;
					threadId = pthread_create(&(treeNode_(tempPos).thread), NULL, &createBddThread, subTree);
					if (threadId != 0){ exit(3); }
				}
				else{ break; }
			}
		}
		currStringPos++; currS = boolFunc[currStringPos];
		oldTreePos = currTreePos;
	}

	if (j == 0){ //incorrect request
		for (i = 1; i < 1024; i++){
			if (treeNode_(i).thread != NULL){
				pthread_join(treeNode_(i).thread, NULL);
			}
		}
		for (i = 1; i < 1024; i++){
			//if (treeNode_(i).thread != NULL){
			//	CloseHandle(treeNode_(i).thread);
			//}
			if (treeNode_(i).bdd != NULL){
				freeBdd((Bdd*)treeNode_(i).bdd);
				treeNode_(i).bdd = NULL;
			}
		}
	}
	else{
		i = treeNode_(0).right;
		subTree = (BTreeSubContainer*)malloc(sizeof(BTreeSubContainer));
		subTree->tree = tree;
		subTree->index = i;
		threadId = pthread_create(&(treeNode_(tempPos).thread), NULL, &createBddThread, subTree);
		pthread_join(treeNode_(i).thread, NULL);
		//CloseHandle(treeNode_(i).thread);
	}

	bdd = (Bdd*)treeNode_(treeNode_(0).right).bdd;
	pthread_mutex_destroy(&(tree->mutex));
	free(tree->node); free(tree);
	return bdd;
}
#endif
