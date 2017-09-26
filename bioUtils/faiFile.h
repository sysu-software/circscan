#ifndef FAI_HEAD_H
#define FAI_HEAD_H

struct faidxInfo {
  int lineLen;
  int lineBlen;
  long len;
  long offset;
};

typedef struct faidxInfo faidx;

typedef map<string, faidx *> faidxMap;

long readFai(FILE *fp, faidxMap &faiHash);

void freeFaiList(faidxMap &fai);

char *faidxFetchSeq(FILE *gfp, const faidx *fai, int start, int end, char strand);

#endif /* End FAI_HEAD_H */
