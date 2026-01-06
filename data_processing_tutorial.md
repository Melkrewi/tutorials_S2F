## data processing
```
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
bamCoverage \
  -b sample.bam \
  -o sample.bw \
  --normalizeUsing CPM \
  --binSize 10

#RNA-seq
bamCoverage -b sample.bam -o sample.bw --normalizeUsing CPM

#paired-end
bamCoverage -b sample.bam -o sample.bw --normalizeUsing CPM --centerReads

#strand-specific RNA-seq
bamCoverage -b sample.bam -o sample.plus.bw  --filterRNAstrand forward
bamCoverage -b sample.bam -o sample.minus.bw --filterRNAstrand reverse

# ChIP-seq / ATAC-seq (often used)
bamCoverage \
  -b sample.bam \
  -o sample.bw \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 2913022398 \
  --extendReads
```
