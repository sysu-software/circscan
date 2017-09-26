#ifndef BED_HEAD_H
#define BED_HEAD_H

struct CBedInfo
{
        struct CBedInfo *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        double score; // normalized read number score
        char strand; // strand
};

typedef struct CBedInfo CBed;

struct CBed6Info
{
        struct CBed6Info *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        char *name;   /* Name of item */
        double score; // normalized read number score
        char strand; // strand
};

typedef struct CBed6Info CBed6;

struct Bed12Info
{
        struct CBed12Info *next; // for list
        char *chrom;    // chromosome name
        int chromStart; // chromosome start
        int chromEnd; // chromosome end
        char *name;     /* Name of item */
        double score; // normalized read number score
        char strand; // strand
        int thickStart; /* Start of where display should be thick (start codon for genes) */
        int thickEnd;   /* End of where display should be thick (stop codon for genes) */
        int itemRgb;    /* RGB 8 bits each */
        int blockCount; /* Number of blocks. */
        int *blockSizes;     /* Comma separated list of block sizes.  */
        int *chromStarts;    /* Start positions inside chromosome.  Relative to chromStart*/
};

typedef struct Bed12Info CBed12;

struct MateBedInfo
{
        struct MateBedInfo *next; // for list
        char *chrom1;    // chromosome name
        int chromStart1; // chromosome start
        int chromEnd1; // chromosome end
        double score1; // normalized read number
        char strand1; // strand
        char *chrom2;    // chromosome name
        int chromStart2; // chromosome start
        int chromEnd2; // chromosome end
        double score2; // normalized read number
        char strand2; // strand
};

typedef struct MateBedInfo MateBed;

void freeBedItem(CBed *bed);

void freeBed6Item(CBed6 *bed6);

void freeBed12Item(CBed12 *bed12);

void freeBedMapList(map<string, CBed *> &bedHash);
/*free bed list */

void freeBedList(CBed *bedList);
/*free bed list */

void sortBed(CBed **pList, int count, int (*compare )(const void *item1,  const void *item2));

void copyBed(CBed *tBed, CBed *oBed);

bool compareBed( const CBedInfo *oBed, const CBedInfo *tBed);

int compareBedScore(const void *a, const void *b);

CBed12 *parseBed12Line(char *line);

CBed6 *parseBed6Line(char *line);

#endif /* End BED_HEAD_H */