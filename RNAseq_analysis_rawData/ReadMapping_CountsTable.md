## Read Mapping and Count table creation for coral larvae against the coral host genome _Acropora digtifera_

The following steps were done for processing the larval RNAseq reads.

1.	Use fastQC v0.11.2 to check the quality of the .fastq files.

```
~/tools/fastqc 26_1_ATCACG_L004_R1_001.fastq
```

2. Remove adapters using cutadapt v3.1.

```
cutadapt -j 4 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 20 -o trimmed1.fastq -p trimmed2.fastq *_R1_001.fastq.gz *_R2_001.fastq.gz
```

3. Make index of the genome.

```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /raid1/home/zoo/kitchens/Adig_genome/bowtie2 --genomeFastaFiles /raid1/home/zoo/kitchens/Adig_genome/bowtie2/adi_v0.9.scaffold.fa --sjdbGTFfile /raid1/home/zoo/kitchens/A.dig_annotate/new_adig.gtf
```

4. Map reads to the genome.

```
STAR --runThreadN 4 --genomeDir /raid1/home/zoo/kitchens/Adig_genome/bowtie2 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix ./26.1_trimmed_ --outBAMsortingBinsN 100 --readFilesIn trimmed1.fastq trimmed2.fastq
```

5. Create counts table with FeatureCounts.
```
featureCounts -p -s 0 -T 8 -a /home/zoo/kitchens/A.dig_annotate/new_adig.gtf -o fc_gene_count_STAR.cnt /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/1rep/26.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/2rep/26.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/3rep/26.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/4rep/26.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/5rep/26.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/6rep/26.6_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/1rep/26.1_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/2rep/26.2_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/3rep/26.3_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/4rep/26.4_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/5rep/26.5_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym26/6rep/26.6_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/1rep/32.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/2rep/32.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/3rep/32.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/4rep/32.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/5rep/32.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/32/6rep/32.6_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/1rep/32.1_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/2rep/32.2_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/3rep/32.3_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/4rep/32.4_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/5rep/32.5_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/sym32/6rep/32.6_sym_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/1rep/72_26.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/2rep/72_26.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/3rep/72_26.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/4rep/72_26.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/5rep/72_26.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/26/6rep/72_26.6_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/1rep/72_sym26.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/2rep/72_sym26.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/3rep/72_sym26.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/4rep/72_sym26.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/5rep/72_sym26.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym26/6rep/72_sym26.6_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/1rep/72_32.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/2rep/72_32.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/3rep/72_32.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/4rep/72_32.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/5rep/72_32.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/32/6rep/72_32.6_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/1rep/72_sym32.1_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/2rep/72_sym32.2_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/3rep/72_sym32.3_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/4rep/72_sym32.4_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/5rep/72_sym32.5_trimmed_Aligned.sortedByCoord.out.bam /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour72/sym32/6rep/72_sym32.6_trimmed_Aligned.sortedByCoord.out.bam
```

6.	Move to R and run deseq2.R script on the counts tables.

## Read mapping symbiotic larval samples against the symbiont, _Symbiodinium tridacnidorum_ (CCMP 2465, clade A3) genome assembly
We also mapped the reads from the larvae exposed to symbionts to the draft genome assembly of _S. tridacnidorum_ to determine if they accounted for some of the reads not mapping to the coral host genome.

1.	Index the _S. tridacnidorum_ genome (clade A3) with STAR.

```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /raid1/home/zoo/kitchens/SymA3/Symbiodinium_tridacnidorum --genomeFastaFiles /raid1/home/zoo/kitchens/SymA3/Symbiodinium_tridacnidorum/Symbiodinium_tridacnidorum.genome.fa --sjdbGTFfile /raid1/home/zoo/kitchens/SymA3/Symbiodinium_tridacnidorum/Symbiodinium_tridacnidorum.gtf
```

2.	The trimmed.fastq file from each symbiotic sample (26 and 32 sym) was run with STAR to align the reads to the genome. 

```
STAR --runThreadN 4 --genomeDir /raid1/home/zoo/kitchens/SymA3 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix ./26.1_A3_trimmed_ --outBAMsortingBinsN 100 --readFilesIn trimmed1.fastq trimmed2.fastq
```

## Phyloflash on all read sets
Due to the lower mapping rates to the host and symbiont genomes, we looked for evidence of other contamination in the raw reads using phyloflash. This tool extracts 16S and 18S ribosomal sequences from metatranscriptomic datasets and compares them against a curated SILVA database. 

Each sample was run as follows:
```
phyloFlash.pl -lib 1d_26_1 -read1 /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/1rep/26_1_ATCACG_L004_R1_001.fastq -read2 /nfs0/IB/Weis_Lab/kitchens/Adigitifera/2014/hour24/26/1rep/26_1_ATCACG_L004_R2_001.fastq -almosteverything -CPUs 2 -dbhome /nfs0/IB/Weis_Lab/kitchens/phyloflash/138
```
