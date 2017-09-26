#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "bioUtils.h"

char **splitWhitespace(char *s, int *entryNum)
/*split the line with whitesapce */
{
int i = 0;
int j = 0;
int len   = strlen(s);
int count = 0;
char **retArray;
char *word;
word = NULL;
retArray = NULL;
while(i < len) {
        while(i < len && isspace(s[i])) ++i;
        j = i;
        while(i < len && !isspace(s[i])) ++i;
        if (j < i) {
                retArray = (char **)safeRealloc(retArray, sizeof(char *)*(count+1));
                word = (char *)safeMalloc(sizeof(char)*(i-j+1));
                strncpy(word, s+j, i-j);
                word[i-j] = '\0';
                retArray[count] = word;
                word = NULL;
                count++;
                }
        }
*entryNum = count;
return retArray;
}

void freeWords(char **s, int entryNum){
int i = 0;
for (i = 0; i < entryNum; i++) {
        safeFree(s[i]);
        }
        safeFree(s);
}

void *safeMalloc(size_t nBytes)
/* safety for malloc */
{
void *p;
p = (void *)calloc(1, nBytes);
if (p == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(EXIT_FAILURE);
        }

return p;
}

void *safeZeroedMalloc(size_t nBytes)
/* safety for malloc. The memory
is initialized to zero. */
{
void *p;
p = (void *)calloc(1, nBytes);
if (p == NULL) {
        fprintf(stderr, "Out of Memory\n");
        exit(EXIT_FAILURE);
        }
memset(p, 0, nBytes);
return p;
}

void *safeRealloc(void *block, size_t nBytes)
/* safety for realloc */
{
void *p;
p = (void *)realloc(block, nBytes);
if (p == NULL) {
        safeFree(block);
        fprintf(stderr, "Out of Memory\n");
        exit(EXIT_FAILURE);
        }
return p;
}

void safeFree(void *block)
/* safety for free() */
{
free(block);
block = NULL;
}

char **splitString(char *s, char *delim, int *entryNum)
/*split the line with whitesapce */
{
int i = 0;
int j = 0;
int len   = strlen(s);
int count = 0;
char **retArray;
char *word;
retArray = NULL;
while(i < len) {
        while(i < len && isDelim(delim, s[i])) ++i;
        j = i;
        while(i < len && !isDelim(delim, s[i])) ++i;
        if (j < i) {
                retArray = (char **)safeRealloc(retArray, (sizeof(char *))*(count+1));
                word = (char *)safeMalloc(sizeof(char)*(i-j+1));
                strncpy(word, s+j, i-j);
                word[i-j] = '\0';
                retArray[count] = word;
                word = NULL;
                count++;
                }
        }
*entryNum = count;
return retArray;
}

boolean isDelim(char *delim, char c)
/* */
{
char *cp;
if ((cp = strchr(delim, c)) == NULL){
        return FALSE;
        }
else {
        return TRUE;
        }
}

char *strClone(char *s)
{
char *p;
p = (char *)safeMalloc(strlen(s)+1);
strcpy(p, s);
return p;
}

char *getLine(FILE *fp)
/* get whole line
 * modified from Vienna RNAfold utils.c 2008/10/4 11:08:01
*/
{
char s[LINE_LEN];
char *line;
char *cp;
boolean done = FALSE;
line = NULL;

do {
    if (fgets(s, LINE_LEN, fp)==NULL)
        break; /* EOF */
    /* for unix OS */
    cp = strchr(s, '\n');
      if (cp != NULL) {
        *cp   = '\0';
        done  = TRUE;
     }
   /* for window OS */
    cp = strchr(s, '\r');
      if (cp != NULL) {
      *cp  = '\0';
      done = TRUE;
     }
    if (line==NULL)
      line = (char *)safeMalloc(strlen(s)+1); /* don't use malloc, for we will using strcat function, so we must initilized the line with '0' */
    else
      line = (char *)safeRealloc(line, strlen(s)+strlen(line)+1);
    strcat(line, s);
  } while(!done);

  return line;
}

char *skipStartWhitespace(char *s)
{
while(isspace(*s)) s++;
        return s;
}

boolean startStr(char *s, char *t)
/* s start with t */
{
        if (!strncmp(s, t, strlen(t))) {
                return TRUE;
                }
        else {
                return FALSE;
                }
}

int overlapLength(int qStart, int qEnd, int tStart, int tEnd)
{
        int maxStart = 0;
        int minEnd = 0;

        maxStart = MAX(qStart, tStart);
        minEnd   = MIN(qEnd, tEnd);
        return (minEnd-maxStart);
}

void toUpperStr(char *string) {
	/* to string upper character */
	long length = strlen(string);
	
	int i;
	for(i=0; i<length; i++) {
		if (isalpha(string[i])){
		string[i] = toupper(string[i]);
	}
		}
	}
  
/*************reverseComp***************/
void reverseComp(char *bytes)
{
/* Reverse complement string */	
long length = strlen(bytes);
reverseBytes(bytes, length);    
complement(bytes);
}

/*************************************/
void complement(char *bytes) {
	   /* complement sequence */
long length = strlen(bytes);
	int i;
for(i=0; i<length; i++) {
	/* bytes[i] = toupper(bytes[i]); */
	switch(bytes[i]){
    case 'A':  bytes[i] = 'T'; break;
    case 'U':  bytes[i] = 'A'; break;
    case 'T':  bytes[i] = 'A'; break;
    case 'C':  bytes[i] = 'G'; break;
    case 'G':  bytes[i] = 'C'; break;
    /* for lower character */
    case 'a':  bytes[i] = 't'; break;
    case 'u':  bytes[i] = 'a'; break;
    case 't':  bytes[i] = 'a'; break;
    case 'c':  bytes[i] = 'g'; break;
    case 'g':  bytes[i] = 'c'; break;    
	 }	 
 }
}

/* Reverse the order of the bytes. */
void reverseBytes(char *bytes, long length)
{
long halfLen = (length>>1);
char *end = bytes+length;
char c;
while (--halfLen >= 0)
    {
    c = *bytes;
    *bytes++ = *--end;
    *end = c;
    }
}

int skipChrom(char *haystack)
{
  const char *needles[3] = {"random", "hap", "gl"};
  int i = 0;
  int returnVal = 0;
  for(i = 0; i < 3; i++)
  {
    //fprintf(stderr, "%s\n", needles[i]);
    if (strstr(haystack, needles[i]))
    {
      return 1;
    }
  }
  return 0;
}

int compIntDescend (const void *a , const void *b)
{
  int *aa = (int * ) a;
  int *bb = (int * )b;
  if( * aa > * bb)return -1;
  if( * aa == * bb) return 0;
  if( * aa < *bb) return 1;
  return 0;
}

int compIntAscend (const void *a , const void *b)
{
  int *aa = (int * )a;
  int *bb = (int * )b;
  if( * aa > * bb)return 1;
  if( * aa == * bb) return 0;
  if( * aa < *bb) return -1;
  return 0;
}

int compDoubleDescend (const void *a , const void *b)
{

  double *aa = (double *)a;
  double *bb = (double *)b;
  if( * aa > * bb)return -1;
  if( * aa == * bb) return 0;
  if( * aa < *bb) return 1;
  return 0;
}

int compDoubleAscend (const void *a , const void *b)
{
  double *aa = (double *)a;
  double *bb = (double *)b;
  if( * aa > * bb)return 1;
  if( * aa == * bb) return 0;
  if( * aa < *bb) return -1;
  return 0;
}

void convertToUpperStr(char *seq)
{
  int k = 0;
  for (k = 0; k < strlen(seq); k++)
  {
    if (seq[k] == 't')
    {
      seq[k] = 'T';
    }
    if (seq[k] == 'a')
    {
      seq[k] = 'A';
    }
    if (seq[k] == 'c')
    {
      seq[k] = 'C';
    }
    if (seq[k] == 'g')
    {
      seq[k] = 'G';
    }
    if (seq[k] == 'u')
    {
      seq[k] = 'U';
    }
  }
}
