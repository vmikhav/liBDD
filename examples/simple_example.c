
#include <stdio.h>
#include "bddcore.h"


int main(){
	Bdd *bdd1, *bdd2, *bdd3;
	char buff[2048];
	gets_s(buff, 2048);
	//example: a1&x1|(a3^1)
	bdd1 = createBdd(buff);
	gets_s(buff, 2048);
	//example: (x1>y)|~z&y
	bdd2 = createBdd(buff);

	bdd3 = synthesizeBdd(bdd1, bdd2, 7, 1, 1);// bdd1|bdd2

	fprintBdd(stdout, bdd1, 0);
	printf("\n");
	fprintBdd(stdout, bdd2, 0);
	printf("\n");
	fprintBdd(stdout, bdd3, 0);

	freeBdd(bdd1); freeBdd(bdd2); freeBdd(bdd3);
	getchar();
}
