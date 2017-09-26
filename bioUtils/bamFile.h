#ifndef BAM_HEAD_H
#define BAM_HEAD_H

struct CSamInfo
{
	struct CSamInfo *next; // for list
	char *chrom;	// chromosome name
	int chromStart; // chromosome start
	int chromEnd; // chromosome end
	int readLen; // read length
	char *readName; // read name
	char *readSeq; // read sequence
	char *mdz; // aligned mdz
	char *cigar; // aligned cigar
	double readNum; // normalized read number
	char strand; // strand
};
typedef struct CSamInfo CSam;

typedef vector<CBed *> bedVector;
typedef map<string, bedVector> chromBedMap;
typedef vector<CSam *> samVector;
typedef map<string, samVector> chromSamMap;

string BuildCigarString(const vector<CigarOp> &cigar);

string PrintTag(const BamAlignment &bam, const string &tag);

void openBamFile(char *bamFile, BamReader &reader);

void openBamNoIdxFile(char *bamFile, BamReader &reader);

int readBamToLocusNum(BamReader &reader, map<string, int> &readLocus);

int readSamToLocusNum(FILE *fp, map<string, int> &readLocus);

int getMDZ(char **fields, int num);

int getEndPos(char *cigar);

double readBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum);

double readBamToSamMap(BamReader &reader, chromSamMap &samHash, map<string, int> &readLocus, int keepDup, int maxLocusNum);

double readPEBamToBedMap(BamReader &reader, chromBedMap &bedHash, map<string, int> &readLocus, int keepDup, int maxLocusNum, int maxInsertLen, int readLen);

int getRegionCounts(BamReader &reader, char *chrom, int chromStart, int chromEnd, char strand);

int filterLowQualities(int qualityBase, int minScore, int minPercent, const char *qualities);

void freeChromBedMap(chromBedMap &bedHash);

void freeBedVector(bedVector &bedList);

void freeSamVector(samVector &samList);

void freeSamItem(CSam *sam);

void freeChromSamMap(chromSamMap &samHash);

bool compareSam(const CSam *oSam, const CSam *tSam);

int compareSamReadNum(const void *a, const void *b);

void copySam(CSam *tSam, CSam *oSam);

#endif /* End BAM_HEAD_H */
