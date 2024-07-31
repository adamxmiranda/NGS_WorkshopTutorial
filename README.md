# NGS Workshop Tutorial

A Resource for the Guo Lab and future NGS analysis

## Types of NGS Sequencing

There are multiple types of NGS assay that each reveal different information. Understanding each assay you wish to analyze is important for choosing the appropriate pipeline of analysis

* ***RNA-seq***
  * for evaluation of gene transcript levels in a cell population
* ***ATAC-seq***
  * for evaluation of accessible chromatin regions
* ***ChIP-seq***
  * Evaluation of the presence and location of DNA bound proteins
* ***HiC***
  * Evaluation of chromatin-chromatin interactions through proximity ligation

## Important File Types to know

General knowledge of the different file types can help you to understand what information is at your disposal and what step in the analytical process you are at.

### FASTQ

This file type is the file type that most sequencers will provide following a sequencing run. These files contain four lines of information for each sequenced read fragment.

general format:

``` markdown
line 1: Read/fragment ID
line 2: The called sequence base pairs (A/G/C/T)
line 3: + and additional identifiers
line 4: quality scores (per base)
```

### SAM/BAM

This file type appends the genomic coordinates to the mapped reads of a fastq file.

### BED

This file type shows a range of genomic coordinates

general format of a bed file:

| Chromosome | Start | End  | Other metadata...
|----------|------------|------------| ------------|
| Chr1   | bp # | bp #|...|
| Chr1    | bp # | bp #|...|
| Chr1    | bp # | bp #|...|


### BEDPE

This file type shows a paired range of coordinates. Usually created for/from HiC experiments to show interacting loci.

| Chromosome | Start | End  | Chromosome | Start | End  |Other metadata...|
|----------|------------|------------| ------------|------------|------------| ------------|
| Chr1   | bp # | bp #| Chr1   | bp # | bp #|...|
| Chr1    | bp # | bp #| Chr1   | bp # | bp #|...|
| Chr1    | bp # | bp #| Chr1   | bp # | bp #|...|


### Counts Matrix

This file contains read counts on a per gene or transcript basis for each sample. in an RNA-seq type assay.

General Format of Counts matrix:

|  | Sample 1 | Sample 2   | Sample 3   | ...|
|----------|------------|------------|------------|------------|
| Gene A   | # of reads | # of reads | # of reads | # of reads |
| Gene B   | # of reads | # of reads | # of reads | # of reads |
| Gene C   | # of reads | # of reads | # of reads | # of reads |
| ...   | # of reads | # of reads | # of reads | # of reads |

### Bigwig

This file type represents the coverage of the genome continuously. Great for the visualization of sequencing data as a continuous histogram.
