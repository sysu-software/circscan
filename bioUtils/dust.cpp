/********************************************************
* the modified code of dust algorithm from NCBI
* dust is A program for filtering low complexity regions 
* from nucleic acid sequences.
* Dust was written by R. Tatusov and D.J. Lipman (unpublished).
 *******************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#include "dust.h"
#include "bioUtils.h"

int dust(int len, char *s, int window, int level, int word)
{
	int i, j, l, from, to, a, b, v, window2;
	int dustNum = 0;
	window2 = window / 2;
	from = 0;
	to = -1;
	for (i=0; i < len; i += window2) {
		from -= window2;
		to -= window2;
		l = (len > i+window) ? window : len-i;
		v = wo(l, s+i, &a, &b, word);
		for (j = from; j <= to; j++) {
			if (isalpha(s[i+j])){
				s[i+j]  = 'N';
				dustNum++;
			}	
		}
		if (v > level) {
			for (j = a; j <= b && j < window2; j++) {
				if (isalpha(s[i+j])){
					s[i+j] = 'N';
					dustNum++;
				}
			}
			from = j;
			to = b;
		} else {
			from = 0;
			to = -1;
		}
	}
	return dustNum;
}


REGION *dust_segs(int len, char *s, int window, int level, int word)
{
	int i, j, l, from, to, a, b, v, nreg, window2;
	REGION *reg = (REGION *)safeMalloc(sizeof(REGION) * MAXREG);;
	nreg = 0;
	from = 0;
	to = -1;
	window2 = window / 2;

	for (i=0; i < len; i += window2) {
		l = (len > i+window) ? window : len-i;
		v = wo(l, s+i, &a, &b, word);
		if (v > level) {
			if (nreg > 0 && reg[nreg-1].to+1 >= a+i) {
				reg[nreg-1].to = b + i;
			} else if (nreg < MAXREG-1) {
				reg[nreg].from = a + i;
				reg[nreg].to = b + i;
				nreg++;
			}
			if (b < window2) {
				i += b - window2;
			}
		}
	}
	reg[nreg].from = 0;
	reg[nreg].to = -1;
	return reg;
}

void wo1(int len, char *s, int ivv, int word)
{
	int i, ii, j, v, t, n, n1, sum;
	int mv, iv, jv;
	static int counts[32*32*32];
	static int iis[32*32*32];
	int js, nis;

	n = 32 * 32 * 32;
	n1 = n - 1;
	nis = 0;
	i = 0;
	ii = 0;
	sum = 0;
	v = 0;
	for (j=0; j < len; j++, s++) {
		ii <<= 5;
		if (isalpha(*s)) {
			if (islower(*s)) {
				ii |= *s - 'a';
			} else {
				ii |= *s - 'A';
			}
		} else {
			i = 0;
			continue;
		}
		ii &= n1;
		i++;
		if (i >= word) {
			for (js=0; js < nis && iis[js] != ii; js++) ;
			if (js == nis) {
				iis[nis] = ii;
				counts[ii] = 0;
				nis++;
			}
			if ((t = counts[ii]) > 0) {
				sum += t;
				v = 10 * sum / j;
				if (mv < v) {
					mv = v;
					iv = ivv;
					jv = j;
				}
			}
			counts[ii]++;
		}
	}
}

int wo(int len, char *s, int *beg, int *end, int word)
{
	int i, j, l1, v;
	int mv, iv, jv;

	l1 = len - word + 1;
	if (l1 < 0) {
		*beg = 0;
		*end = len - 1;
		return 0;
	}
	mv = 0;
	iv = 0;
	jv = 0;
	for (i=0; i < l1; i++) {
		wo1(len-i, s+i, i, word);
	}
	*beg = iv;
	*end = iv + jv;
	return mv;
}
