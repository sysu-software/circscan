#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
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
#include "faiFile.h"
#include "bamFile.h"

template<typename T>
string NumberToString(T Number)
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

string BuildCigarString(const vector<CigarOp> &cigar)
{

  stringstream cigarString;
  for (size_t i = 0; i < cigar.size(); ++i)
  {
    switch (cigar[i].Type)
    {
    case ('M') :
    case ('I') :
    case ('D') :
    case ('N') :
    case ('S') :
    case ('H') :
    case ('P') :
      cigarString << cigar[i].Length << cigar[i].Type;
    }
  }
  return cigarString.str();
}

string PrintTag(const BamAlignment &bam, const string &tag)
{
  uint32_t uTagValue;
  int32_t sTagValue;
  ostringstream value;
  if (bam.GetTag(tag, uTagValue))
    value << uTagValue;
  else if (bam.GetTag(tag, sTagValue))
    value << sTagValue;
  else
  {
    cerr << "The requested tag ("
       << tag
       << ") was not found in the BAM file.  Exiting\n";
    exit(1);
  }
  return value.str();
}

void openBamFile(char *bamFile, BamReader &reader)
{
  if (!reader.Open(bamFile))
  {
    cerr << "Failed to open BAM file " << bamFile << endl;
    exit(1);
  }
  else
  {
    reader.LocateIndex();
  }
  if ( !reader.HasIndex())
  {
    cerr << "Failed to open BAI file " << bamFile << endl;
    //exit(1);
  }
}

void openBamNoIdxFile(char *bamFile, BamReader &reader)
{
  if (!reader.Open(bamFile))
  {
    cerr << "Failed to open BAM file " << bamFile << endl;
    exit(1);
  }
}

int readBamToLocusNum(BamReader &reader, map<string, int> &readLocus)
{
  int i = 0;
  reader.Rewind();
  
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName = bam.Name;
      readLocus[readName] += 1;
      i++;
    }
  }
  return i;
}


int readSamToLocusNum(FILE *fp, map<string, int> &readLocus)
{
  char *strLine;
  char *tmpLine;
  char delims[] = "-";
  int headTag = 0;
  int fieldNum = 0;
  char **fields   = NULL;
  char **infos = NULL;
  int infoNum = 0;
  int i = 0;
  int flag = 0;
  double totalNum = 1;
  if (feof(fp))
  {
    return 0;
  }
  //chr1  147646  147749  U6-related  1000  -
  while(strLine = getLine(fp))
  {
    tmpLine  = strLine;
    tmpLine  = skipStartWhitespace(tmpLine);
    fieldNum = 0;
    fields   = NULL;
    if (tmpLine[0] == '@')
    {
      safeFree(strLine);
      continue;
    }

    fields = splitWhitespace(tmpLine, &fieldNum);
    if(fieldNum >= 10 && strcmp(fields[2], "*") != 0)
    {
      string readName(fields[0]);
      if (readLocus.find(readName) == readLocus.end())
      {
        readLocus[readName] = 1;
      }
      else
      {
        readLocus[readName] += 1;
      }
      i++;
    }
    freeWords(fields, fieldNum);
    safeFree(strLine);
  }
  return i;
}

double readBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CBed *bedPtr    = NULL;
  map<string, int> dupMap;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum) continue; // remove multiple locus
      string chrom  = refs.at(bam.RefID).RefName;
      char strand = '+';
      if (bam.IsReverseStrand() == true) strand = '-';
      string dupPos = chrom + NumberToString(bam.Position)+ NumberToString(bam.GetEndPosition(false, false))+strand;
      if (!keepDup)
      {
        if (dupMap.find(dupPos) != dupMap.end())
        {
          dupMap[dupPos] += 1;
          continue;
        }
        else{
          dupMap[dupPos]  = 1;
        }
      }
      bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      bedPtr->strand     = strand;
      //if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      bedPtr->score    = 1 / (double)lociNum;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          bedPtr->score = atof(infos[1]) / (double)lociNum;
      }
      totalNum += bedPtr->score;
      freeWords(infos, infoNum); // free memory
      bedPtr->next = NULL;
      bedHash[chrom].push_back(bedPtr);
      i++;
    }
  }
  return totalNum;
}

int filterLowQualities(int qualityBase, int minScore, int minPercent, const char *qualities)
{
  int i = 0;
  int count = 0;
  int seqLen = strlen(qualities);
  for(i=0; i<seqLen; i++){
    int val = qualities[i]-qualityBase-minScore;
    //fprintf(stderr, "%c %d %d\n", qualities[i], qualities[i], val);
    if (val > 0) count++;
  }
  if (count/(double)seqLen*100>=minPercent) return 1;
  return 0;
}

double readPEBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, 
       int keepDup, int maxLocusNum, int maxInsertLen, int readLen)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  int matePosition = 0;
  int insertSize   = 0;
  int read1Len     = 0;
  int read2Len     = 0;
  int maxDiffLen   = 5;
  map<string, int> pairHash;
  CBed *bedPtr     = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum * 2) continue; // remove multiple locus
      /*const char* qualities = bam.Qualities.c_str();
      if (filterLowQualities(33, 20, 90, qualities) == 0){
        continue;
      }*/
      bedPtr             = (CBed *)safeMalloc(sizeof(CBed));
      string chrom       = refs.at(bam.RefID).RefName;
      bedPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      bedPtr->chromStart = bam.Position;
      bedPtr->chromEnd   = bam.GetEndPosition(false, false);
      read1Len           = bedPtr->chromEnd - bedPtr->chromStart;
      bedPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) bedPtr->strand = '-';
      if (bam.IsSecondMate() == true)
      {
        if (bedPtr->strand == '-')
        {
          bedPtr->strand = '+';
        }
        else
        {
          bedPtr->strand = '-';
        }
      }
      matePosition = bam.MatePosition;
      insertSize = bam.InsertSize;
      if (insertSize != 0 && bam.RefID == bam.MateRefID && abs(insertSize)<maxInsertLen)
      {
        if (read1Len<readLen && abs(bedPtr->chromStart-matePosition)>maxDiffLen){
          freeBedItem(bedPtr);
          continue;
        }
        string pairString = readName + NumberToString(bedPtr->chromStart);
        if (pairHash.find(pairString) == pairHash.end())
        {
          string readString = readName + NumberToString(matePosition);
          pairHash[readString] = read1Len;
          freeBedItem(bedPtr);
          continue;
        }
        else{
          read2Len = pairHash[pairString];
        }
        if (insertSize > 0)
        {
          bedPtr->chromEnd = bedPtr->chromStart + insertSize;
        }
        else
        {
          bedPtr->chromStart = matePosition;
        }
      }// if insert size
      else
      {
        freeBedItem(bedPtr);
        continue;
      }
      if ((bedPtr->chromEnd-bedPtr->chromStart)> readLen 
           && (read1Len<readLen-maxDiffLen || read2Len<readLen-maxDiffLen)) // discard the degradome reads
      {
        freeBedItem(bedPtr);
        continue;
      }
      bedPtr->score = 1.0 / (double)lociNum;
      if (lociNum > 1) bedPtr->score *= 2;
      totalNum += bedPtr->score;
      bedPtr->next = NULL;
      bedHash[chrom].push_back(bedPtr);
      i++;
    }
  }
  return totalNum;
}

double readBamToSamMap(BamReader &reader, chromSamMap &samHash, map<string, int> &readLocus, int keepDup, int maxLocusNum)
{
  int i           = 0;
  char delims[]   = "-";
  int fieldNum    = 0;
  char **infos    = NULL;
  int infoNum     = 0;
  double totalNum = 0;
  CSam *samPtr    = NULL;
  // get header & reference information
  reader.Rewind();
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  // rip through the BAM file and convert each mapped entry to BED
  BamAlignment bam;
  while (reader.GetNextAlignment(bam))
  {
    if (bam.IsMapped() == true)
    {
      string readName    = bam.Name;
      int lociNum        = 1;
      if (readLocus.find(readName) != readLocus.end())
      {
        lociNum = readLocus[readName];
      }
      if (lociNum > maxLocusNum) continue; // remove multiple locus
      samPtr             = (CSam *)safeMalloc(sizeof(CSam));
      string chrom       = refs.at(bam.RefID).RefName;
      samPtr->chrom      = strClone(const_cast<char *>(chrom.c_str()));
      samPtr->chromStart = bam.Position;
      samPtr->chromEnd   = bam.GetEndPosition(false, false);
      samPtr->readName   = strClone(const_cast<char*>(readName.c_str()));
      samPtr->strand     = '+';
      if (bam.IsReverseStrand() == true) samPtr->strand = '-';
      string cigar    = BuildCigarString(bam.CigarData);
      samPtr->cigar   = strClone(const_cast<char*>(cigar.c_str()));
      string readSeq  = bam.QueryBases;
      samPtr->readSeq = strClone(const_cast<char*>(readSeq.c_str()));
      string mdz;
      if (!bam.GetTag("MD", mdz)) {
        fprintf(stderr, "The requested tag MD is not exist in bam file\n");
        exit(1);
      }
      samPtr->mdz = strClone(const_cast<char*>(mdz.c_str()));
      samPtr->readNum    = 1 / (double)lociNum;
      infos = splitString(const_cast<char *>(readName.c_str()), delims, &infoNum);
      if (infoNum == 2 && keepDup)
      {
        if (isdigit(infos[1][0]))
          samPtr->readNum = atof(infos[1]) / (double)lociNum;
      }
      totalNum += samPtr->readNum;
      freeWords(infos, infoNum); // free memory
      samPtr->next = NULL;
      samHash[chrom].push_back(samPtr);
      i++;
    }
  }
  return totalNum;
}

int getRegionCounts(BamReader &reader, char *chrom, int chromStart, int chromEnd, char strand)
{
  int regionReadNum = 0;
  reader.Rewind();
  int id = reader.GetReferenceID(chrom);
  BamRegion region(id, chromStart, id, chromEnd);
  // rip through the BAM file
  if ( (id != -1) && (reader.SetRegion(region)) )
  {
    BamAlignment bam;
    while (reader.GetNextAlignment(bam))
    {
      if (bam.IsMapped() == true)
      {
        char bamStrand   = '+';
        if (bam.IsReverseStrand() == true) bamStrand = '-';
        if (strand != bamStrand) continue;
        int regStart = bam.Position;
        int regEnd = bam.GetEndPosition(false, false);
        if (overlapLength(chromStart, chromEnd, regStart, regEnd)>0)
        {
            regionReadNum += 1;
        }
      } // is mapped
    } //while bam
  } // if id AND
  return regionReadNum;
}


int getMDZ(char **fields, int num)
{
  int i = 0;
  int mdzIdx = 0;
  const char *mdz = "MD:Z:";
  for (i = 0; i < num; i++)
  {
    if (strncmp(fields[i], mdz, 5) == 0)
    {
      mdzIdx = i;
      break;
    }
  }
  return mdzIdx;
}

int getEndPos(char *cigar)
{
  const char *samcigar = "MIDNSHP=X";
  int start = 0;
  int end   = 0;
  int readLen   = 0;
  int i = 0;
  char c;
  for (i = 0; i < strlen(cigar); i++)
  {
    c = cigar[i];
    if (!isdigit(c))
    {
      cigar[i] = '\0';
      int len = atoi(cigar + start);
      switch ( c )
      {
      // increase end position on CIGAR chars [DMXN=]
      case 'D' :
      case 'M' :
      case 'X' :
      case 'N' :
      case '=' :
        readLen += len;
        break;
      }
      cigar[i] = c;
      start = i + 1;
    }
  }
  return readLen;
}

void freeBedVector(bedVector &bedList)
{
    for (bedVector::iterator vecItr = bedList.begin(); vecItr != bedList.end(); vecItr++)
    {
      CBed *bed = *vecItr;
      freeBedItem(bed);
    }
    bedList.clear();
}

void freeSamVector(samVector &samList)
{
    for (samVector::iterator vecItr = samList.begin(); vecItr != samList.end(); vecItr++)
    {
      CSam *sam = *vecItr;
      freeSamItem(sam);
    }
    samList.clear();
}


void freeChromBedMap(chromBedMap &bedHash) /*free bed map */
{
  for (chromBedMap::iterator mapItr = bedHash.begin(); mapItr != bedHash.end(); mapItr++) {
    bedVector bedList = mapItr->second;
    freeBedVector(bedList);
  }
  bedHash.clear();
}

void freeChromSamMap(chromSamMap &samHash) /*free sam map */
{
  for (chromSamMap::iterator mapItr = samHash.begin(); mapItr != samHash.end(); mapItr++) {
    samVector samList = mapItr->second;
    freeSamVector(samList);
  }
  samHash.clear();
}

void freeSamItem(CSam *sam)
{
  safeFree(sam->chrom);
  safeFree(sam->readName);
  safeFree(sam->readSeq);
  safeFree(sam->mdz);
  safeFree(sam->cigar);
  safeFree(sam);
}

void copySam(CSam *tSam, CSam *oSam)
// copy oBed to tBed
{
  tSam->chrom      = strClone(oSam->chrom);
  tSam->chromStart = oSam->chromStart;
  tSam->chromEnd   = oSam->chromEnd;
  tSam->readLen    = oSam->readLen;
  tSam->readLen    = oSam->readLen;
  tSam->readNum    = oSam->readNum;
  tSam->readSeq    = strClone(oSam->readSeq);
  tSam->readName   = strClone(oSam->readName);
  tSam->mdz        = strClone(oSam->mdz);
  tSam->cigar      = strClone(oSam->cigar);
}

bool compareSam(const CSam *oSam, const CSam *tSam)
{
  return (oSam->chromStart<tSam->chromStart);
}

int compareSamReadNum(const void *a, const void *b)
{
  CSam *oSam = *(CSam **)a;
  CSam *tSam = *(CSam **)b;
  if ((tSam->readNum - oSam->readNum) > 0)
    return 1;
  else if ((tSam->readNum - oSam->readNum) < 0)
    return -1;
  else
    return 0;
}
