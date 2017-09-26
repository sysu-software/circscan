/*****************************************************************************************
 *	circScan - An algorithm for discovering circRNAs from short-read CLIP-seq data
 *
 *	Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>
 * 
 *	School of Life Sciences, Sun Yat-Sen University
 *	
 *  Create date: 09/18/2015
 *  
 ****************************************************************************************/

#ifndef circScan_HEAD_H
#define circScan_HEAD_H

#ifndef MIN
#define MIN(x, y) (x)<(y)?(x):(y)
#endif

#ifndef MAX
#define MAX(x, y) (x)>(y)?(x):(y)
#endif

#define GAP_OPEN -2.0
#define GAP_CONT -2.0
#define MATCH 1.0
#define MISMATCH -2.0

#define H 3

struct parameterInfo
{
	int verbose;
	int minSeqLen;
	int maxMismatch;
	int seedLen;
	int parclip;
	int maxNum;
	int minScore;
	double errorRate;
	double fdr;
};

typedef struct parameterInfo parameters;

struct acceptorIndexInfo
{
	char *name;
	char strand;
	int  chromPos;
	int  seqStart;
	int  exonIdx;
	int  misPos;
	int  t2cNum;
};
typedef struct acceptorIndexInfo acceptorIdxInfo;

struct donorIndexInfo
{
	int  chromPos;
	int  seqStart;
	int  exonIdx;
	int  misPos;
	int  t2cNum;
};
typedef struct donorIndexInfo donorIdxInfo;

struct bedSeqInfo
{
	char *chrom;
	char *seq;
};
typedef struct bedSeqInfo seqInfo;

struct circJunctionInfo
{
	char *chrom;
	int chromStart;
	int chromEnd;
	char *name;
	int score;
	char strand;
	int mismatchNum;
	int junctionPos;
	int t2cNum;
	char *readName;
	char *readSeq;
};

typedef struct circJunctionInfo junctionInfo;

typedef map<string, seqInfo *> seqMap;
typedef vector<acceptorIdxInfo *> acceptorVector;
typedef map<string, acceptorVector> acceptorIdxMap;
typedef map<string, donorIdxInfo *> donorMap;
typedef map<string, donorMap> donorIdxMap;
typedef vector<junctionInfo *> junctionVector;

void scanCRI(parameters *paraInfo, FILE *outfp, FILE *genomefp, FILE *faifp, FILE *readfp, FILE *bedfp, char *bamFile);
void getSeqIndex(parameters *paraInfo, FILE *genomefp, faidx *fai, CBed12 *bed12, seqMap &seqHash,
				 acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash);
void searchCircRNAs(parameters *paraInfo, FILE *outfp, BamReader &reader, FILE *fp, seqMap &seqHash,
					acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash);
int searchJunctions(parameters *paraInfo, FILE *outfp, BamReader &reader, char *seq, char *readName,
					seqMap &seqHash, acceptorIdxMap &acceptorIdxHash, donorIdxMap &donorIdxHash,
					junctionVector &jVector);
int extendMatchSeq(parameters *paraInfo, char *readSeq, char *readName, seqInfo *sif, acceptorIdxInfo *aidx,
				   donorIdxInfo *didx, int junctionPos, junctionVector &jVector);
void outputHeader(FILE *outfp);
void outputJunctionVector(parameters *paraInfo, FILE *outfp, BamReader &reader, junctionVector &jList, double fdr);
int maskDustNum(char *seq, int seqLen);

int skipNoSplit(char *dSite, char *aSite, char strand);

void copyAcceptorInfo(acceptorIdxInfo *org, acceptorIdxInfo *obj);

void copyDonorInfo(donorIdxInfo *org, donorIdxInfo *obj);

void freeJunction(junctionInfo *jInfo);

void freeJunctionVector(junctionVector &jList);

void freeSeqMap(seqMap &seqHash);

void freeSeqInfo(seqInfo *sif);

void freeAcceptorIdxMap(acceptorIdxMap &acceptorIdxHash);

void freeAcceptorVector(acceptorVector &aList);

void freeAcceptorIdxItem(acceptorIdxInfo *al);

void freeDonorIdxMap(donorIdxMap &donorIdxHash);

void freeDonorMap(donorMap &dHash);

void freeDonorIdxItem(donorIdxInfo *didx);

int randNum(int len);

int seqShuffle(char *s1, char *s2);

#endif /* End circScan_HEAD_H */
