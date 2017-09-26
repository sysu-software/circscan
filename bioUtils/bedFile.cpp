#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"
#include "bedFile.h"

bool compareBed( const CBedInfo *oBed, const CBedInfo *tBed)
{
  return (oBed->chromStart < tBed->chromStart);
}

int compareBedScore(const void *a, const void *b)
{
  CBed *oBed = *(struct CBedInfo **)a;
  CBed *tBed = *(struct CBedInfo **)b;
  if ((tBed->score - oBed->score) > 0)
    return 1;
  else if ((tBed->score - oBed->score) < 0)
    return -1;
  else
    return 0;
}

void freeBedList(CBed *bedList)
/*free bed list */
{
  CBed *m = NULL;
  CBed *p = NULL;
  m = bedList;
  while(m != NULL)
  {
    p = m;
    m = p->next;
    freeBedItem(p);
  }
}

void freeBedItem(CBed *bed)
{
  safeFree(bed->chrom);
  safeFree(bed);
}

void freeBed6Item(CBed6 *bed6)
{
  safeFree(bed6->chrom);
  safeFree(bed6->name);
  safeFree(bed6);
}

void freeBed12Item(CBed12 *bed12)
{
  safeFree(bed12->chrom);
  safeFree(bed12->name);
  safeFree(bed12->blockSizes);
  safeFree(bed12->chromStarts);
  safeFree(bed12);
}


void freeBedMapList(map<string, CBed *> &bedHash) /*free bed list */
{
  for (map<string, CBed *>::iterator curr = bedHash.begin(); curr != bedHash.end();
      curr++)
  {
    CBed *p = curr->second;
    freeBedList(p);
  }
  bedHash.clear();
}


void copyBed(CBed *tBed, CBed *oBed)
// copy oBed to tBed
{
  tBed->chrom      = strClone(oBed->chrom);
  tBed->chromStart = oBed->chromStart;
  tBed->chromEnd   = oBed->chromEnd;
  tBed->score    = oBed->score;
  tBed->strand     = oBed->strand;
}

void sortBed(CBed **pList, int count, int (*compare )(const void *item1,  const void *item2))
// sort bed, I have done major modified from jksrc common.c slSort, thanks Jim
{
  //fprintf(stderr, "sort bed array...\n");
  CBed *list = *pList;
  if (count > 1)
  {
    CBed *el;
    CBed **array;
    int i;
    array = (CBed **)safeMalloc(count * sizeof(CBed *));
    for (el = list, i = 0; el != NULL; el = el->next, i++)
      array[i] = el;
    qsort(array, count, sizeof(array[0]), compare);
    list = array[0];
    for (i = 1; i < count; ++i)
    {
      list->next = array[i];
      list = array[i];
    }
    list->next = NULL;
    *pList = array[0];
    safeFree(array);
  }
  //fprintf(stderr, "sort end...\n");
}

CBed12 *parseBed12Line(char *line)
{
  int i = 0;
  int fieldNum = 0;
  char **fields = NULL;
  int tmpFieldNum = 0;
  char **tmpFields = NULL;
  char delim[] = ",";
  CBed12 *bed = NULL;
  fields = splitWhitespace(line, &fieldNum);
  if (fieldNum != 12)
  {
    freeWords(fields, fieldNum);
    fprintf(stderr, "the format of annotation file must be bed12\n");
    exit(1);
  }
  bed = (CBed12 *)safeMalloc(sizeof(CBed12));
  bed->chrom = strClone(fields[0]);
  bed->chromStart = atoi(fields[1]);
  bed->chromEnd = atoi(fields[2]);
  bed->name = strClone(fields[3]);
  bed->score = atof(fields[4]);
  bed->strand = fields[5][0];
  bed->thickStart = atoi(fields[6]);
  bed->thickEnd = atoi(fields[7]);
  bed->itemRgb = 0;
  bed->blockCount = atoi(fields[9]);
  bed->blockSizes = (int *)safeMalloc(sizeof(int) * bed->blockCount);
  bed->chromStarts = (int *)safeMalloc(sizeof(int) * bed->blockCount);
  tmpFields = splitString(fields[10], delim, &tmpFieldNum);
  for(i = 0; i < bed->blockCount; i++)
  {
    bed->blockSizes[i] = atoi(tmpFields[i]);
  }
  freeWords(tmpFields, tmpFieldNum);
  tmpFields = splitString(fields[11], delim, &tmpFieldNum);
  for(i = 0; i < bed->blockCount; i++)
  {
    bed->chromStarts[i] = atoi(tmpFields[i]);
  }
  freeWords(tmpFields, tmpFieldNum);
  freeWords(fields, fieldNum);

  return bed;
}

CBed6 *parseBed6Line(char *line)
{
  int i = 0;
  int fieldNum = 0;
  char **fields = NULL;
  char delim[] = ",";
  CBed6 *bed = NULL;
  fields = splitWhitespace(line, &fieldNum);
  if (fieldNum < 6)
  {
    freeWords(fields, fieldNum);
    fprintf(stderr, "the format of annotation file must be be6\n");
    exit(1);
  }
  bed = (CBed6 *)safeMalloc(sizeof(CBed6));
  bed->chrom = strClone(fields[0]);
  bed->chromStart = atoi(fields[1]);
  bed->chromEnd = atoi(fields[2]);
  bed->name = strClone(fields[3]);
  bed->score = atof(fields[4]);
  bed->strand = fields[5][0];
  freeWords(fields, fieldNum);
  return bed;
}
