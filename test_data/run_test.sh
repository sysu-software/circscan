#!/usr/bin/env bash

echo ">>> building bowtie2 index..."
bowtie2-build testGenome.fa testGenome
echo ">>> building faidx..."
samtools faidx testGenome.fa
echo ">>> aligning example reads"
bowtie2 -p 16 -t -k 4 --no-unal -D 200 -R 3 -N 0 -L 15 -i S,1,0.5 --score-min=C,-16,0 -f -U testReadFile.fa  -x testGenome --un ./testGenome.unmapped.fa
echo ">>> converting alignments to bam"
samtools view -bS testGenome.bwt2 >testGenome.bam
echo ">>> sorting bam"
samtools sort testGenome.bam testGenome.sorted
echo ">>> building bam index"
samtools index testGenome.sorted.bam
echo ">>> runing the circScan to get circRNAs"
../bin/circScan testGenome.fa testGenome.fa.fai testAnnoFile.bed12 testGenome.sorted.bam testGenome.unmapped.fa >test_circScan_candidates.txt

rm testGenome.bwt2
