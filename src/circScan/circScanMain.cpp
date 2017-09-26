/*****************************************************************************************
 *  circScan - An algorithm for discovering circRNAs from short-read CLIP-seq data
 *
 *  Author : Jian-Hua Yang <yangjh7@mail.sysu.edu.cn>
 * 
 *  School of Life Sciences, Sun Yat-Sen University
 *  
 *  Create date: 09/18/2015
 *  
 ****************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include "BamReader.h"
#include "BamAux.h"
using namespace BamTools;
#include <getopt.h>
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
#include "circScan.h"

void usage(void);
char version[] = "circScan version 0.1 [2015/11/06]";

int main(int argc, char **argv)
{
  int i           = 0;
  int c           = 0;
  char *outfile   = NULL;
  char *bamFile   = NULL;
  FILE *outfp     = NULL;
  FILE *genomefp  = NULL;
  FILE *faifp     = NULL;
  FILE *readfp    = NULL;
  FILE *bedfp     = NULL;
  int showVersion = 0;
  int showHelp    = 0;
  parameters paraInfo;

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVdpa:o:m:x:s:e:i:f:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "par-clip" , no_argument , NULL, 'p' },
    { "outfile" , required_argument , NULL, 'o' },
    { "mismatch" , required_argument , NULL, 'x' },
    { "seed-len" , required_argument , NULL, 's' },
    { "error-rate" , required_argument , NULL, 'e' },
    { "max-num" , required_argument , NULL, 'm' },
    { "min-score" , required_argument , NULL, 'i' },
    { "fdr" , required_argument , NULL, 'f' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose     = 0;
  paraInfo.minSeqLen   = 12;
  paraInfo.maxMismatch = 1;
  paraInfo.seedLen     = 12;
  paraInfo.parclip     = 0;
  paraInfo.errorRate   = 0.05;
  paraInfo.maxNum      = 1;
  paraInfo.minScore    = 20;
  paraInfo.fdr         = 0.05;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'p':
      paraInfo.parclip = 1;
      break;
    case 'o':
      outfile  = optarg;
      break;
    case 'm':
      paraInfo.maxNum = atoi(optarg);
      break;
    case 'x':
      paraInfo.maxMismatch = atoi(optarg);
      break;
    case 's':
      paraInfo.seedLen = atoi(optarg);
      break;
    case 'e':
      paraInfo.errorRate = atof(optarg);
      break;
    case 'i':
      paraInfo.minScore = atoi(optarg);
      break;
    case 'f':
      paraInfo.fdr = atof(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (paraInfo.minSeqLen < 10)
  {
    fprintf(stderr, "-m option:  minimum sequence length (for each anchor) must be large than 10 nt\n");
    usage();
  }
  if (paraInfo.seedLen < 10)
  {
    fprintf(stderr, "-s option:  minimum seed length must be large than 10 nt\n");
    usage();
  }
  if (paraInfo.maxMismatch > 2)
  {
    fprintf(stderr, "-x option: maximum mismatch number is 2\n");
    usage();
  }

  if (argc == optind || argc - 5 != optind) usage();
  genomefp = (FILE *) fopen(argv[optind], "r");
  if (genomefp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open your genome file: %s\n", argv[optind]);
    usage();
  }
  faifp = (FILE *) fopen(argv[optind + 1], "r");
  if (faifp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open your fai file: %s\n", argv[optind + 1]);
    usage();
  }
  bedfp = (FILE *) fopen(argv[optind + 2], "r");
  if (bedfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open your annotation file: %s\n", argv[optind + 2]);
    usage();
  }
  bamFile = argv[optind + 3];
  readfp = (FILE *) fopen(argv[optind + 4], "r");
  if (readfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open your read file: %s\n", argv[optind + 3]);
    usage();
  }
  if (outfile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outfile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outfile);
      usage();
    }
  }

  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  fprintf(stderr, "# circScan program start\n");
  scanCRI(&paraInfo, outfp, genomefp, faifp, readfp, bedfp, bamFile);
  fprintf(stderr, "# circScan program end\n");
  fclose(faifp);
  fclose(genomefp);
  fclose(outfp);

  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  circScan [options] <genome file> <genome fai> <annotation file, bed12> <mapped alignments> <unmapped read>\n\
File format for annotation must be UCSC bed12\n\
Mapped alignments must be sort BAM format\n\
[options]\n\
-v/--verbose                :verbose information\n\
-V/--version                : rseScan version\n\
-o/--outfile <string>       : output file\n\
-p/--parclip                : par-clip datasets\n\
-f/--fdr                    : false discovery rate(<0.05)\n\
-m/--max-num <int>          : maximum locations for each read [default=1]\n\
-s/--seed-len <int>         : the length of seed for searching [default>=10]\n\
-x/--mismatch <int>         : maximum mismatch number for each read [default<=1]\n\
-e/--error-rate <double>    : maximum error rate [default<=0.05], =mismatch/readLen\n\
-i/--min-score <int>        : minimum score [default>=20], =matches - mismatches\n\
");

  exit(1);
}
