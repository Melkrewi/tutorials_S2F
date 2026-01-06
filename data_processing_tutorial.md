## data processing

In a slurm script: `fastqc.sh`

Content of the script:
```ruby
#load modules
module load java
module load fastqc

#command
srun fastqc *.fastq.gz
```

Look at html files 

# TRIMMOMATIC Single End (SE) 

### Forward reads (*_R1_001.fastq.gz)

Use the `HEADCROP` option trimm the first X bp of the reads (here the first first 30 bp of the reads for example)

In a slurm script: `trimmomatic_SE_forward.sh`

Content of the script:

```ruby
#load modules
module load java

#command
srun java -jar /PATH_TO/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 40 SAMPLEID_R1_001.fastq.gz SAMPLEID_R1_001_SE.fastq.gz HEADCROP:30
```

### Reverse reads (*R2_001.fastq.gz)

Use the `HEADCROP` option trimm the first X bp of the reads (here the first first 8 bp of the reads for example)

In a slurm script: `trimmomatic_SE_reverse.sh`

Content of the script:
```ruby
#load modules
module load java

#command
srun java -jar /PATH_TO/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 40 SAMPLEID_R2_001.fastq.gz SAMPLEID_R2_001_SE.fastq.gz HEADCROP:8
```

# TRIMMOMATIC Paired End (PE) 

Slurm script : `trimmomatic_PE.sh`

Content of the slurm script:
```ruby
#load modules
module load java

#command
srun java -jar /PATH_TO/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 40 SAMPLEID_R1_001_SE.fastq.gz SAMPLEID_R2_001_SE.fastq.gz SAMPLEID_R1_001_PE_paired.fastq.gz SAMPLEID_R1_001_PE_unpaired.fastq.gz SAMPLEID_R2_001_PE_paired.fastq.gz SAMPLEID_R2_001_PE_unpaired.fastq.gz ILLUMINACLIP:/PATH_TO/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:22 MINLEN:36 LEADING:3 TRAILING:3 
```

Comments: Look up these parameters in the [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)
- ILLUMINACLIP:/PATH_TO/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10
- SLIDINGWINDOW:4:22
- MINLEN:36
- LEADING:3
- TRAILING:3

# Read trimming with Trimgalore

## You will need to:
* load trimgalore, cutadapt and fastqc
* run Trimgalore

## Add to your SLURM script:

```
#load modules here
module load trimgalore/0.5.0 
module load cutadapt 
module load fastqc/0.11.7  

#run commands here
srun trim_galore --fastqc --paired -q 20 --length 35 reads_1.fastq reads_2.fastq
```

[Main change relative to default: keep only reads longer than 35 post-trimming (20 seems short for SNP calling, as it may lead to misalignments):

## Short explanation of the options used:
--fastqc: Run FastQC in the default mode on the FastQ file once trimming is complete.

--paired: This option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files. To pass the validation test, both sequences of a sequence pair are required to have a certain minimum length which is governed by the option --length (see above). If only one read passes this length threshold the other read can be rescued (see option --retain_unpaired). Using this option lets you discard too short read pairs without disturbing the sequence-by-sequence order of FastQ files which is required by many aligners.

--length <INT>: Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. Default: 20 bp.

-q/--quality <INT>: Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal). Default Phred score: 20.]



# Read mapping:
### list of tools
* Bowtie2
* BWA
* Hisat2
* STAR
* Segemehl
* Subread

## Bowtie2:
```
module load bowtie2/2.4.5
srun bowtie2-build fasta.fna sicusGenome 
srun bowtie2 -x sicusGenome -U 317397_S19_R1_001.fastq --local -p 50 -S female8.sam

```
## BWA MEM:
```
module load bwa
module load samtools
bwa index GENOME_ASSEMBLY.fasta

for i in *_1.fastq
do
    prefix=$(basename $i _1.fastq)
    srun bwa mem -M -t 60 GENOME_ASSEMBLY.fasta ${prefix}_1.fastq ${prefix}_2.fastq | samtools view -@30 -F 4 -b | samtools sort -@30 -T individual > ${prefix}.bam
    srun samtools index -@30 ${prefix}.bam
done
```
## STAR
```
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /path/to/STAR_index \
  --genomeFastaFiles genome.fa \
  --sjdbGTFfile annotation.gtf \
  --sjdbOverhang 149
STAR \
  --runThreadN 8 \
  --genomeDir /path/to/STAR_index \
  --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix sample_ \
  --outSAMtype BAM SortedByCoordinate
```
## Segemehl
```
module load conda
conda activate segemehl
segemehl.x -x reference.idx -d GENOME_ASSEMBLY.fasta -t 60 #index genome
segemehl.x -i reference.idx -d GENOME_ASSEMBLY.fasta -q READS_1_paired.fastq -p READS_5_2_paired.fastq -t 60 > alignment.sam #align reads
module load samtools
samtools view -@ 30 -bS -o alignment.bam alignment.sam #convert sam to bam
samtools sort -@ 30 -o alignment.sorted.bam alignment.bam #sort bam
samtools index alignment.sorted.bam #index bam
rm alignment.sam
rm alignment.bam
samtools view -@ 25 -F0x100 -b -o alignment.sorted_no_secondary.bam alignment.sorted.bam   #remove secondary alignments
samtools index -@25 alignment.sorted_no_secondary.bam #index bam
```
## Hisat2
```
module load hisat2
module load stringtie
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_asexuality_redo_2/2.Transcriptomic_short_reads_alignments/
mkdir hisat2
hisat2-build GENOME_ASSEMBLY.fasta genome_index
for i in *_1.fastq.gz
do
   prefix=$(basename $i _1.fastq.gz)
   hisat2 --phred33 -p 50 --novel-splicesite-outfile hisat2/${prefix}_splicesite.txt -S hisat2/${prefix}_accepted_hits.sam -x genome_index -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz --rna-strandness RF --max-intronlen 50000
   samtools view -@ 25 -bS -o hisat2/${prefix}_accepted_hits.bam hisat2/${prefix}_accepted_hits.sam
   samtools sort -@ 25 -o hisat2/${prefix}_accepted_hits.sorted.bam hisat2/${prefix}_accepted_hits.bam
   samtools view -@ 25 -F0x100 -bS -h hisat2/${prefix}_accepted_hits.sorted.bam -o hisat2/${prefix}_aln_no_secondary.bam
   samtools sort -@ 25 -o hisat2/${prefix}_aln_no_secondary.sorted.bam hisat2/${prefix}_aln_no_secondary.bam
   samtools index hisat2/${prefix}_aln_no_secondary.sorted.bam
done
```
## Subread
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_Nauplii_project_round_2/21.subread_package/subread-2.0.2-Linux-x86_64/bin/:$PATH

subread-buildindex -o genome_index GENOME_ASSEMBLY.fasta

subread-align -T 20 -t 1 -d 50 -D 600 -i genome_index -r READS_1.fastq.gz -R READS_2.fastq.gz -o subread_results.bam
```


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
