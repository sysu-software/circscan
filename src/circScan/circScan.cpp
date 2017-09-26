/********************************************************
 * circScan: search for circRNA-RBP interactions
 * from CLIP-seq datasets
 * jianhua yang yangjh7@mail.sysu.edu.cn
 * $ 2015/11/6 Sun Yat-sen University
 *******************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
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
#include "dust.h"
#include "circScan.h"

template<typename T>
string NumberToString(T Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}

void scanCRI(parameters *paraInfo, FILE *outfp, FILE *genomefp, FILE *faifp,
			 FILE *readfp, FILE *bedfp, char *bamFile)
/* scan circRNA-RBP interactions */
{
	BamReader reader;
	char *line = NULL;
	CBed12 *bed12 = NULL;
	seqMap seqHash;
	acceptorIdxMap acceptorIdxHash;
	donorIdxMap donorIdxHash;
	faidxMap faiHash;
	if(paraInfo->verbose) fprintf(stderr, "read bam file\n");
	openBamFile(bamFile, reader);
	if(paraInfo->verbose) fprintf(stderr, "read genome fai file\n");
	readFai(faifp, faiHash);
	if (paraInfo->verbose) fprintf(stderr, "index annotations\n");
	outputHeader(outfp);
	while (line = getLine(bedfp))
	{
		if (feof(bedfp) || line == NULL)
		{
			safeFree(line);
			break;
		}
		bed12 = parseBed12Line(line); // get bed information
		safeFree(line); // free memory

		string chrom(bed12->chrom);
		if (faiHash.find(chrom) == faiHash.end())
		{
			fprintf(stderr, "can't not find the chromosome %s, skip it.", bed12->chrom);
			freeBed12Item(bed12);
			continue;
		}
		faidx *fai = faiHash[chrom];
		getSeqIndex(paraInfo, genomefp, fai, bed12, seqHash, acceptorIdxHash, donorIdxHash);
		freeBed12Item(bed12);
	}// loop-while
	searchCircRNAs(paraInfo, outfp, reader, readfp, seqHash, acceptorIdxHash, donorIdxHash);
	// free all memory
	freeSeqMap(seqHash);
	freeAcceptorIdxMap(acceptorIdxHash);
	freeDonorIdxMap(donorIdxHash);
	freeFaiList(faiHash);
}

void outputHeader(FILE *outfp)
{
	fprintf(outfp, "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tmismatchNum\tt2cNum\tacceptReadNum\tdonorReadNum\tFDR\tjunctionPosInRead\treadName\treadSeq\n");
	fflush(outfp);
}

void searchCircRNAs(parameters *paraInfo, FILE *outfp, BamReader &reader, FILE *fp, seqMap &seqHash,
					acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash)
{
	int fieldNum      = 0;
	char **fields     = NULL;
	char *line        = NULL;
	char *seq         = NULL;
	char *qualityName = NULL;
	char *quality     = NULL;
	char readName[255];
	int readLen       = 0;
	// for dust parameters
	int level         = 20;
	int word          = 3;
	int window        = 64;

	if (paraInfo->verbose) fprintf(stderr, "search circRNAs\n");
	while (line = getLine(fp))
	{
		if (feof(fp) || line == NULL)
		{
			safeFree(line);
			break;
		}
		if (line[0] != '@' && line[0] != '>')
		{
			fprintf(stderr, "error read format: %c\n", line[0]);
			safeFree(line);
			continue;
		}
		if (line[0] == '@')
		{
			fields      = splitWhitespace(line + 1, &fieldNum);
			seq         = getLine(fp);
			qualityName = getLine(fp);
			quality     = getLine(fp);
			strcpy(readName, fields[0]);
			freeWords(fields, fieldNum);
			safeFree(qualityName);
			safeFree(quality);
		}
		else if (line[0] == '>')
		{
			fields = splitWhitespace(line + 1, &fieldNum);
			seq    = getLine(fp);
			strcpy(readName, fields[0]);
			freeWords(fields, fieldNum);
		}
		readLen = strlen(seq);
		if (dust(readLen, seq, window, level, word) < 1 && readLen >= paraInfo->seedLen * 2)
		{
			junctionVector jVector;
			int jTag = searchJunctions(paraInfo, outfp, reader, seq, readName, seqHash, acceptorIdxHash, donorIdxHash, jVector);
			int shuffleNum = 1000;
			int passNum = 0;
			if (jTag)
			{
				int i = 0;
				int seqLen = strlen(seq);
				char *shuffleSeq = (char *)safeMalloc(sizeof(char) * (seqLen + 1));
				strcpy(shuffleSeq, seq);
				for(i = 0; i < shuffleNum; i++)
				{
					junctionVector sVector;
					seqShuffle(shuffleSeq, seq);
					//fprintf(stderr, "%s\n", shuffleSeq);
					if(searchJunctions(paraInfo, outfp, reader, shuffleSeq, readName, seqHash, acceptorIdxHash, donorIdxHash, sVector) > 0)
					{
						passNum++;
					}
					freeJunctionVector(sVector);
				}
				safeFree(shuffleSeq);
				double fdr = (double)passNum / (double)shuffleNum;
				if (fdr < paraInfo->fdr) outputJunctionVector(paraInfo, outfp, reader, jVector, fdr);
			}
			freeJunctionVector(jVector);
		}
		safeFree(seq);
		safeFree(line);
	} // gofile while loops
}

int maskDustNum(char *seq, int seqLen)
{
	int i = 0;
	for (i = 0; i < seqLen; i++)
	{
		if (seq[i] == 'N') return 0;
	}
	return 1;
}

int searchJunctions(parameters *paraInfo, FILE *outfp, BamReader &reader, char *seq, char *readName,
					seqMap &seqHash, acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash,
					junctionVector &jVector)
{
	int i = 0;
	int j = 0;
	int seqLen = strlen(seq);
	int startPos = paraInfo->seedLen;
	int endPos = seqLen - paraInfo->seedLen;
	int mapNum = 0;
	int spliceSignal = 0;
	map<string, int> aMap; // remove duplications
	map<string, int> dMap; // remove duplications

	//if (paraInfo->verbose) fprintf(stderr, "search for junctions...\n");
	char *acceptorSeq = (char *)safeMalloc(sizeof(char) * (paraInfo->seedLen + 1));
	char *donorSeq = (char *)safeMalloc(sizeof(char) * (paraInfo->seedLen + 1));
	for (i = startPos; i < endPos; i++)
	{
		strncpy(acceptorSeq, seq + i, paraInfo->seedLen);
		strncpy(donorSeq, seq + i - paraInfo->seedLen, paraInfo->seedLen);
		acceptorSeq[paraInfo->seedLen] = '\0';
		donorSeq[paraInfo->seedLen] = '\0';
		string aiseq(acceptorSeq);
		string diseq(donorSeq);
		if (acceptorIdxHash.find(aiseq) != acceptorIdxHash.end())
		{
			acceptorVector avector = acceptorIdxHash[aiseq];
			for (acceptorVector::iterator vecItr = avector.begin(); vecItr != avector.end(); vecItr++)
			{
				acceptorIdxInfo *aidx = *vecItr;
				string geneName(aidx->name);
				seqInfo *sif = seqHash[geneName];
				string aKey = sif->chrom + NumberToString(aidx->chromPos);
				if (aMap.find(aKey) != aMap.end()) continue; // remove duplications
				donorMap dmap = donorIdxHash[geneName];
				if (dmap.find(diseq) != dmap.end())
				{
					donorIdxInfo *didx = dmap[diseq];
					string dKey = sif->chrom + NumberToString(didx->chromPos);
					if (dMap.find(dKey) != dMap.end()) continue; // remove duplications
					//if (paraInfo->verbose) fprintf(stderr, "%s\t%d\t%d\t%c\t%s\t%s\t%s\n", aidx->name, aidx->chromPos, didx->chromPos, aidx->strand, acceptorSeq, donorSeq, sif->seq);
					if (aidx->exonIdx <= didx->exonIdx)
					{
						int retVal = extendMatchSeq(paraInfo, seq, readName, sif, aidx, didx, i, jVector);
						if (retVal == 2) spliceSignal = 1;
						if (retVal && spliceSignal == 0)
						{
							aMap[aKey] = aidx->chromPos; // for removing duplications
							dMap[dKey] = didx->chromPos; // for removing duplications
						}
					} // exonidx
				} // dmap
			} // iterator
		} // sequences
	}// for startPos-endPos
	mapNum = jVector.size();
	safeFree(donorSeq);
	safeFree(acceptorSeq);
	if (mapNum > 0 && mapNum <= paraInfo->maxNum) return 1;
	return 0;
}

void outputJunctionVector(parameters *paraInfo, FILE *outfp, BamReader &reader, junctionVector &jList, double fdr)
{
	int extendLen = 100;
	for (junctionVector::iterator vecItr = jList.begin(); vecItr != jList.end(); vecItr++)
	{
		junctionInfo *jInfo = *vecItr;
		int regStart = jInfo->chromStart;
		int regEnd = regStart + extendLen;
		int aReadNum = getRegionCounts(reader, jInfo->chrom, regStart, regEnd, jInfo->strand);
		regEnd = jInfo->chromEnd;
		regStart = regEnd - extendLen;
		if (regStart < 0) regStart = 0;
		int dReadNum = getRegionCounts(reader, jInfo->chrom, regStart, regEnd, jInfo->strand);
		fprintf(outfp, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\t%d\t%.5f\t%d\t%s\t%s\n",
				jInfo->chrom, jInfo->chromStart, jInfo->chromEnd, jInfo->name,
				jInfo->score, jInfo->strand, jInfo->mismatchNum, jInfo->t2cNum,
				aReadNum, dReadNum, fdr,
				jInfo->junctionPos, jInfo->readName, jInfo->readSeq);
		fflush(outfp);
	}
}

int extendMatchSeq(parameters *paraInfo, char *readSeq, char *readName, seqInfo *sif, acceptorIdxInfo *aidx,
				   donorIdxInfo *didx, int junctionPos, junctionVector &jVector)
{
	//if (paraInfo->verbose) fprintf(stderr, "extending sequences... %d %d %d\n", aidx->seqStart, didx->seqStart, junctionPos);
	int i = 0;
	int j = 0;
	int sucess = 1;
	int failure = 0;
	int mismatchNum = 0;
	int t2cNum = 0;
	char *geneSeq = sif->seq;
	int readLen = strlen(readSeq);
	int geneLen = strlen(geneSeq);
	if (geneLen - aidx->seqStart < readLen - junctionPos) return failure;
	if (didx->seqStart < junctionPos) return failure;
	j = aidx->seqStart + paraInfo->seedLen;
	for (i = junctionPos + paraInfo->seedLen; i < readLen && i < geneLen; i++)
	{
		if (readSeq[i] != geneSeq[j])
		{
			if (paraInfo->parclip && ((readSeq[i] == 'C' || readSeq[i] == 'c') && (geneSeq[j] == 'T' || geneSeq[j] == 't'))) // for par-clip experiments
			{
				t2cNum += 1;
			}
			else
			{
				mismatchNum += 1;
			}
		}
		if (mismatchNum > paraInfo->maxMismatch) return failure;
		j++;
	}
	j = didx->seqStart - paraInfo->seedLen;
	for (i = junctionPos - paraInfo->seedLen; i >= 0; i--)
	{
		if (readSeq[i] != geneSeq[j])
		{
			if (paraInfo->parclip && ((readSeq[i] == 'C' || readSeq[i] == 'c') && (geneSeq[j] == 'T' || geneSeq[j] == 't'))) // for par-clip experiments
			{
				t2cNum += 1;
			}
			else
			{
				mismatchNum += 1;
			}
		}
		if (mismatchNum > paraInfo->maxMismatch) return failure;
		j--;
	}
	if (mismatchNum / readLen > paraInfo->errorRate) return failure;

	if (readLen - mismatchNum * 2 < paraInfo->minScore) return failure;

	if (aidx->strand == '+' && aidx->chromPos > didx->chromPos)
	{
		sucess = 2;
		return sucess;
	}
	if (aidx->strand == '-' && aidx->chromPos < didx->chromPos)
	{
		sucess = 2;
		return sucess;
	}
	int chromStart = aidx->chromPos;
	int chromEnd = didx->chromPos;
	if (aidx->chromPos > didx->chromPos)
	{
		chromStart = didx->chromPos;
		chromEnd = aidx->chromPos;
	}
	junctionInfo *jInfo = (junctionInfo *)safeMalloc(sizeof(junctionInfo));
	jInfo->chrom = strClone(sif->chrom);
	jInfo->chromStart = chromStart;
	jInfo->chromEnd = chromEnd;
	jInfo->name = strClone(aidx->name);
	jInfo->score = readLen - mismatchNum * 2;
	jInfo->strand = aidx->strand;
	jInfo->junctionPos = junctionPos;
	jInfo->mismatchNum = mismatchNum;
	jInfo->t2cNum = t2cNum + aidx->t2cNum + didx->t2cNum;
	jInfo->readName = strClone(readName);
	jInfo->readSeq = strClone(readSeq);

	jVector.push_back(jInfo);
	if (paraInfo->verbose) fprintf(stderr, "#%s\t%d\t%d\t%s\t%d\t%c\t%s\t%s\n", sif->chrom, chromStart, chromEnd, aidx->name,
									   readLen - mismatchNum, aidx->strand, readName, readSeq);
	return sucess;
}

int skipNoSplit(char *dSite, char *aSite, char strand)
{
	int returnVal = 0;
	if (strand == '-')
	{
		char *tmp = dSite;
		dSite = aSite;
		aSite = tmp;
	}
	if (toupper(dSite[0]) != 'G' || toupper(dSite[1]) != 'T' || toupper(aSite[0]) != 'A' || toupper(aSite[1]) != 'G')
	{
		returnVal = 1;
	}
	return returnVal;
}

void getSeqIndex(parameters *paraInfo, FILE *genomefp, faidx *fai, CBed12 *bed12, seqMap &seqHash,
				 acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash)
{
	int i = 0;
	int j = 0;
	int idx = 0;
	int geneLen = 0;
	int acceptorPos = 0;
	int donorPos = 0;
	char *seq;
	char *acceptorSeq = (char *)safeMalloc(sizeof(char) * (paraInfo->seedLen + 1));
	char *donorSeq = (char *)safeMalloc(sizeof(char) * (paraInfo->seedLen + 1));
	donorMap donorHash;
	acceptorIdxInfo *acceptorIdx = NULL;
	donorIdxInfo *donorIdx = NULL;
	for(i = 0; i < bed12->blockCount; i++)
	{
		geneLen += bed12->blockSizes[i];
	}
	seq = (char *)safeMalloc(sizeof(char) * (geneLen + 1));
	for(i = 0; i < bed12->blockCount; i++)
	{
		if (bed12->strand == '+')
		{
			idx = i;
			if(idx > 0) acceptorPos += bed12->blockSizes[idx - 1]; // plus location
		}
		else  // reverse the all items
		{
			idx = bed12->blockCount - i - 1;
			if(idx < bed12->blockCount - 1) acceptorPos += bed12->blockSizes[idx + 1]; // minus location
		}
		donorPos += bed12->blockSizes[idx];
		int chromStart = bed12->chromStart + bed12->chromStarts[idx];
		int chromEnd = chromStart + bed12->blockSizes[idx];
		char *subSeq = faidxFetchSeq(genomefp, fai, chromStart, chromEnd, bed12->strand);
		strcat(seq, subSeq);
		if (bed12->blockSizes[idx] < paraInfo->minSeqLen || bed12->blockSizes[idx] < paraInfo->seedLen)
		{
			safeFree(subSeq);
			continue;
		}
		char *dSite = faidxFetchSeq(genomefp, fai, chromEnd, chromEnd + 2, bed12->strand);
		char *aSite = faidxFetchSeq(genomefp, fai, chromStart - 2, chromStart, bed12->strand);
		if(skipNoSplit(dSite, aSite, bed12->strand))
		{
			safeFree(dSite);
			safeFree(aSite);
			safeFree(subSeq);
			continue;
		}
		acceptorIdx           = (acceptorIdxInfo *)safeMalloc(sizeof(acceptorIdxInfo));
		donorIdx              = (donorIdxInfo *)safeMalloc(sizeof(donorIdxInfo));
		acceptorIdx->name     = strClone(bed12->name);
		acceptorIdx->strand   = bed12->strand;
		acceptorIdx->exonIdx  = i;
		acceptorIdx->seqStart = acceptorPos;
		acceptorIdx->misPos   = -1;
		acceptorIdx->t2cNum   = 0;
		donorIdx->exonIdx     = i;
		donorIdx->seqStart    = donorPos;
		donorIdx->misPos   	  = -1;
		donorIdx->t2cNum      = 0;
		if (bed12->strand == '+')
		{
			acceptorIdx->chromPos = chromStart;
			donorIdx->chromPos = chromEnd;
			strncpy(acceptorSeq, subSeq, paraInfo->seedLen);
			strncpy(donorSeq, subSeq + bed12->blockSizes[idx] - paraInfo->seedLen, paraInfo->seedLen);
		}
		else  // reverse the all items
		{
			acceptorIdx->chromPos = chromEnd;
			donorIdx->chromPos = chromStart;
			strncpy(acceptorSeq, subSeq, paraInfo->seedLen);
			strncpy(donorSeq, subSeq + bed12->blockSizes[idx] - paraInfo->seedLen, paraInfo->seedLen);
		}
		acceptorSeq[paraInfo->seedLen] = '\0';
		donorSeq[paraInfo->seedLen] = '\0';
		string aiseq(acceptorSeq);
		acceptorIdxHash[aiseq].push_back(acceptorIdx);
		string diseq(donorSeq);
		donorHash[diseq] = donorIdx;
		if (paraInfo->parclip)
		{
			for (j = 0; j < paraInfo->seedLen; j++)
			{
				if (acceptorSeq[j] == 'T' || acceptorSeq[j] == 't')
				{
					acceptorIdxInfo *acceptorIdx2  = (acceptorIdxInfo *)safeMalloc(sizeof(acceptorIdxInfo));
					copyAcceptorInfo(acceptorIdx, acceptorIdx2);
					acceptorSeq[j] = 'C';
					acceptorIdx2->misPos = j;
					acceptorIdx2->t2cNum = 1;
					string aiseq(acceptorSeq);
					acceptorIdxHash[aiseq].push_back(acceptorIdx2);
					acceptorSeq[j] = 'T';
				}
				if (donorSeq[j] == 'T' || donorSeq[j] == 't')
				{
					donorIdxInfo *donorIdx2 = (donorIdxInfo *)safeMalloc(sizeof(donorIdxInfo));
					copyDonorInfo(donorIdx, donorIdx2);
					donorSeq[j] = 'C';
					donorIdx2->misPos = j;
					donorIdx2->t2cNum = 1;
					string diseq(donorSeq);
					donorHash[diseq] = donorIdx2;
					donorSeq[j] = 'T';
				}
			}
		}
		// free memory
		safeFree(dSite);
		safeFree(aSite);
		safeFree(subSeq);
	}
	seq[geneLen] = '\0';
	string mapKey(bed12->name);
	seqInfo *sif = (seqInfo *)safeMalloc(sizeof(seqInfo));
	sif->seq = seq;
	sif->chrom = strClone(bed12->chrom);
	seqHash[mapKey] = sif;
	donorIdxHash[mapKey] = donorHash;
	safeFree(donorSeq);
	safeFree(acceptorSeq);
}

void copyAcceptorInfo(acceptorIdxInfo *org, acceptorIdxInfo *obj)
{
	obj->name     = strClone(org->name);
	obj->strand   = org->strand;
	obj->exonIdx  = org->exonIdx;
	obj->seqStart = org->seqStart;
	obj->chromPos = org->chromPos;
}

void copyDonorInfo(donorIdxInfo *org, donorIdxInfo *obj)
{
	obj->exonIdx  = org->exonIdx;
	obj->seqStart = org->seqStart;
	obj->chromPos = org->chromPos;
}

void freeJunction(junctionInfo *jInfo)
{
	safeFree(jInfo->chrom);
	safeFree(jInfo->name);
	safeFree(jInfo->readName);
	safeFree(jInfo->readSeq);
	safeFree(jInfo);
}

void freeJunctionVector(junctionVector &jList)
{
	for (junctionVector::iterator vecItr = jList.begin(); vecItr != jList.end(); vecItr++)
	{
		junctionInfo *junc = *vecItr;
		freeJunction(junc);
	}
	jList.clear();
}

void freeDonorIdxMap(donorIdxMap &donorIdxHash)
{
	for (donorIdxMap::iterator mapItr = donorIdxHash.begin(); mapItr != donorIdxHash.end(); mapItr++)
	{
		donorMap dHash = mapItr->second;
		freeDonorMap(dHash);
	}
	donorIdxHash.clear();
}

void freeDonorMap(donorMap &dHash)
{
	for (donorMap::iterator mapItr = dHash.begin(); mapItr != dHash.end(); mapItr++)
	{
		donorIdxInfo *didx = mapItr->second;
		freeDonorIdxItem(didx);
	}
	dHash.clear();
}

void freeDonorIdxItem(donorIdxInfo *didx)
{
	safeFree(didx);
}

void freeAcceptorIdxMap(acceptorIdxMap &acceptorIdxHash)
{
	for (acceptorIdxMap::iterator mapItr = acceptorIdxHash.begin(); mapItr != acceptorIdxHash.end(); mapItr++)
	{
		acceptorVector avector = mapItr->second;
		freeAcceptorVector(avector);
	}
	acceptorIdxHash.clear();
}


void freeAcceptorVector(acceptorVector &aList)
{
	for (acceptorVector::iterator vecItr = aList.begin(); vecItr != aList.end(); vecItr++)
	{
		acceptorIdxInfo *al = *vecItr;
		freeAcceptorIdxItem(al);
	}
	aList.clear();
}

void freeAcceptorIdxItem(acceptorIdxInfo *al)
{
	safeFree(al->name);
	safeFree(al);
}

void freeSeqMap(seqMap &seqHash)
{
	for (seqMap::iterator mapItr = seqHash.begin(); mapItr != seqHash.end(); mapItr++)
	{
		seqInfo *sif = mapItr->second;
		freeSeqInfo(sif);
	}
	seqHash.clear();
}

void freeSeqInfo(seqInfo *sif)
{
	safeFree(sif->chrom);
	safeFree(sif->seq);
}

int randNum(int len)
{
	if (len <= 0)
	{
		fprintf(stderr, "random seed must be >=1");
		exit(1);
	}
	int randVal = 0;
	srand((unsigned)time(NULL));
	randVal = (int)(rand() % len);
	return randVal;
}

int seqShuffle(char *s1, char *s2)
{
	int  len;
	int  pos;
	char c;

	if (s1 != s2) strcpy(s1, s2);
	for (len = strlen(s1); len > 1; len--)
	{
		pos       = randNum(len);
		//fprintf(stderr, "%d\n", pos);
		c         = s1[pos];
		s1[pos]   = s1[len - 1];
		s1[len - 1] = c;
	}
	return 1;
}

