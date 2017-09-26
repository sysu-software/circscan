/********************************************************
* the modified code of dust algorithm from NCBI
* dust is A program for filtering low complexity regions 
* from nucleic acid sequences.
* Dust was written by R. Tatusov and D.J. Lipman (unpublished).
 *******************************************************/
#ifndef DUST_HEAD_H
#define DUST_HEAD_H

#define MAXREG    1001

typedef struct {
	int from;
	int to;
	int score;
} REGION;

REGION *dust_segs(int len, char *s, int window, int level, int word);

void wo1(int len, char *s, int ivv, int word);

int wo(int len, char *s, int *beg, int *end, int word);

int dust(int len, char *s, int window, int level, int word);


#endif /* End DUST_HEAD_H */