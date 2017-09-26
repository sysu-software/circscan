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
#include "faiFile.h"

long readFai(FILE *fp, faidxMap &faiHash)
{
  char *line = NULL;
  char *tmpLine = NULL;
  int fieldNum = 0;
  char **fields = NULL;
  char delim[] = "\t";
  faidx *fai = NULL;
  long genomeSize = 0;
  while (line = getLine(fp))
  {
    tmpLine = line;
    tmpLine = skipStartWhitespace(tmpLine);
    if (feof(fp) || tmpLine == NULL)
    {
      safeFree(tmpLine);
      break;
    }
    fields = splitString(tmpLine, delim, &fieldNum);
    if (fieldNum < 5)
    {
      continue;
    }
    string chrom(fields[0]);
    fai            = (faidx *) safeMalloc(sizeof(faidx));
    fai->len       = atol(fields[1]);
    fai->offset    = atol(fields[2]);
    fai->lineBlen  = atoi(fields[3]);
    fai->lineLen   = atoi(fields[4]);
    genomeSize     += fai->len;
    faiHash[chrom] = fai;
    freeWords(fields, fieldNum);
    safeFree(line);
  } // gofile while loops
  return genomeSize;
}

void freeFaiList(faidxMap &fai)
/*free fai list */
{
  for (faidxMap::iterator curr = fai.begin(); curr != fai.end();
      curr++)
  {
    safeFree(curr->second);
  }
}

char *faidxFetchSeq(FILE *gfp, const faidx *fai, int start, int end, char strand)
{
  int l;
  char c;
  char *seq = NULL;
  l = 0;
  if (start < 0)
    start = 0;
  if (end < 0)
    end = 0;
  if (start > fai->len)
    start = fai->len;
  if (end > fai->len)
    end = fai->len;
  seq = (char *)malloc(end - start + 1);
  long offset = fai->offset + int(start / fai->lineBlen) * fai->lineLen + start % fai->lineBlen;
  fseek(gfp, offset, SEEK_SET);
  while (fread(&c, 1, 1, gfp) == 1 && l < end - start)
    if (isgraph(c)) seq[l++] = toupper(c);
  seq[l] = '\0';
  if(strand == '-')
  {
    reverseComp(seq);
  }
  return seq;
}
