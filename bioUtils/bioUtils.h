#ifndef UTILS_HEAD_H
#define UTILS_HEAD_H

typedef int boolean;
#define FALSE 0
#define TRUE  1

#define EXIT_FAILURE 1
#define EXIT_SUCCESS 0

#define LINE_LEN 512

#ifndef MIN
#define MIN(x, y) (x)<(y)?(x):(y)
#endif

#ifndef MAX
#define MAX(x, y) (x)>(y)?(x):(y)
#endif

/* Some stuff to support large files in Linux. */
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE 1
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif

extern char **splitWhitespace(char *s, int *entryNum);
/*split the line with whitesapce */

extern void freeWords(char **s, int entryNum);
/* free words */

extern void *safeMalloc(size_t bytes);
/* safety for malloc */

extern void *safeZeroedMalloc(size_t nBytes);

extern void *safeRealloc(void *block, size_t bytes);
/* safety for realloc */

extern void safeFree(void *block);
/* safety for free() */

extern char **splitString(char *s, char *delim, int *entryNum);
/*split the line with whitesapce */

extern boolean isDelim(char *delim, char c);
/* */

extern char *strClone(char *s);
/* clone string */

extern char *getLine(FILE *fp);
/* get whole line
 * modified from Vienna RNAfold utils.c 2008/10/4 11:08:01
*/

extern char *skipStartWhitespace(char *s);
/* skip leading whitespace */

extern boolean startStr(char *s, char *t);
/* s start with t */

extern int overlapLength(int qStart, int qEnd, int tStart, int tEnd);
/* caculate the overlap length */

/* follows function from Jim Kent jksrc/common.h  Thanks Jim */
extern void reverseBytes(char *bytes, long length);
/* reverse bytes */

extern void reverseComp(char *bytes);
/* reverse complement */

extern void complement(char *bytes);
/* complement sequence */

extern void toUpperStr(char *string);

int skipChrom(char *haystack);

int compIntDescend (const void *a ,const void *b);

int compIntAscend (const void *a ,const void *b);

int compDoubleDescend (const void *a ,const void *b);

int compDoubleAscend (const void *a ,const void *b);

void convertToUpperStr(char *seq);

#endif /* End UTILS_HEAD_H */
