## Read Mapping and Count table creation for coral larvae against the coral host genome _Acropora digtifera_

The following steps were done for processing the larval RNAseq reads.

1.	Index the _A. digitifera_ genome with Bowtie2.

```
/local/cluster/bin/bowtie2-build /home/zoo/kitchens/A.dig_annotate/adi_v0.9.scaffold.fa adi_v0.9.scaffold
```

2.	Use fastQC v0.11.2 to check the quality of the .fastq files.

```
~/tools/fastqc 26_1_ATCACG_L004_R1_001.fastq
```

3.	Remove adaptor sequences using cutadapt v1.6 as follows:

```
~/tools/cutadapt -a AGATCGGAAGAGC --minimum-length 20 -o trimmed1.fastq -p trimmed2.fastq 26_1_ATCACG_L004_R1_001.fastq 26_1_ATCACG_L004_R2_001.fastq
```

4.	The trimmed.fastq file from each sample was run with TopHat v2.0.12 to align the reads to the genome. 

tophat –p (number of processors) –G (annotated.gff3 file) –o (output directory) (genome fasta file) (input file names)

Ex.
```
/local/cluster/bin/tophat -p 8 -G /home/zoo/kitchens/A.dig_annotate/new_adig.gff3 -o tophat_26.1_trimmed /home/zoo/kitchens/Adig_genome/bowtie2/adi_v0.9.scaffold trimmed1.fastq trimmed2.fastq
```

5.	The number of mapped reads in the accepted_hits.bam file was tabulated using samtools. These values are reported in Table S1. 
```
samtools flagstat ./tophat_26.1_trimmed/accepted_hits.bam
```

6.	The mapped reads in thee accepted_hits.bam file were converted raw counts using samtools and bedtools as follows. 

Each BAM file was sorted and indexed.
```
#Sort:
/local/cluster/bin/samtools sort ./tophat_26.1_trimmed/accepted_hits.bam sort26.1

# Index:
/local/cluster/bin/samtools index sort26.1.bam
```
Then, read counts were extracted from all sorted-indexed BAM files using all the gene features of the A. digitifera GFF file and combined using bedtools multicov tool.
```
# Make table:
bedtools multicov -bams sort26.1.bam sort26.2.bam sort26.3.bam sort26.4.bam sort26.5.bam sort26.6.bam \
sort26.1sym.bam sort26.2sym.bam sort26.3sym.bam sort26.4_sym.bam sort26.5_sym.bam  sort26.6sym.bam \
sort32.1.bam sort32.2.bam sort32.3.bam sort32.4.bam sort32.5.bam sort32.6.bam \
sort32.1sym.bam  sort32.2sym.bam sort32.3sym.bam sort32.4sym.bam sort32.5sym.bam sort32.6sym.bam \
72sort26.1.bam 72sort26.2.bam 72sort26.3.bam 72sort26.4.bam 72sort26.5.bam 72sort26.6.bam \
72sort26.1sym.bam 72sort26.2sym.bam 72sort26.3sym.bam 72sort26.4sym.bam 72sort26.5sym.bam 72sort26.6sym.bam \
72sort32.1.bam 72sort32.2.bam 72sort32.3.bam 72sort32.4.bam 72sort32.5.bam 72sort32.6.bam \
72sort32.1sym.bam 72sort32.2sym.bam 72sort32.3sym.bam 72sort32.4sym.bam 72sort32.5sym.bam 72sort32.6sym.bam \
-bed /home/zoo/kitchens/A.dig_annotate/new_adig.gff3 > all_counts.gff
```

The read counts for each gene were extracted and cleaned up.
```
# Pull out “gene” only:
grep -P "\tgene\t" all_counts.gff > gene_only.gff

# Truncate table to gene name and count data:
sed 's/^.*Name=//' gene_only.gff > gene_counts.tab
```

7.	Move to R and run deseq2.R script on the counts table.

## Read mapping symbiotic larval samples against the symbiont, _Symbiodinium tridacnidorum_ (CCMP 2465, clade A3) genome assembly
We also mapped the reads from the larvae exposed to symbionts to the draft genome assembly of _S. tridacnidorum_ to determine if they accounted for some of the reads not mapping to the coral host genome.

1.	Index the _S. tridacnidorum_ genome (clade A3) with Bowtie2.

```
/local/cluster/bin/bowtie2-build symA3_140711_k2.final.ov1k.removedup.fasta symA3_140711_k2.final.ov1k.removedup.fasta
```

2.	The trimmed.fastq file from each symbiotic sample (26 and 32 sym) was run with TopHat to align the reads to the genome. There was no annotation file for this genome at the time these samples were processed (2014). 

```
/local/cluster/bin/tophat -p 8 -o SymA3_tophat_26.1sym /nfs0/IB/Weis_Lab/kitchens/symA_genome/symA3_140711_k2.final.ov1k.removedup trimmed1.fastq trimmed2.fastq
```

## Phyloflash on all read sets
Due to the lower mapping rates to the host and symbiont genomes, we looked for evidence of other contamination in the raw reads using phyloflash. This tool extracts 16S and 18S ribosomal sequences from metatranscriptomic datasets and compares them against a curated SILVA database. 

Each sample was run as follows:
```
phyloFlash.pl -lib 1d_26_1 -read1 /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/1rep/26_1_ATCACG_L004_R1_001.fastq -read2 /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/1rep/26_1_ATCACG_L004_R2_001.fastq -almosteverything -CPUs 2 -dbhome /nfs0/IB/Weis_Lab/kitchens/phyloflash/138
```