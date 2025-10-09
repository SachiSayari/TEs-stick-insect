# TEs-stick-insect
Work with NGS raw reads, dnaPipeTE

# DATA FRAME
_Bacillus_ _spp._


|SPECIES	  |TAX. ID	|REPRODUCTION
| -------- | ------- | --------- |
|*B. grandii*	|55088		|facultative p. |
|*B. rossius*	|75184		|obligate b. |
|*B. atticus*	|36825		|obligate p. |


# DATA
RNA-seq reads accession numbers: SRX7034623 to SRX7034670 (Bioproject: PRJNA578804)    
Assembled transcriptomes accession numbers GJDY01000000, GJDZ01000000, GJEA00000000    

# USED COMMANDS
First, we used **prefetch** in a loop to download SRA data for every species into a separate folder. Folders `B_atticus, B_rossius, B_grandii` with accession no's in ```.tsv```.
Example of used command:
```
while read -r acc; do prefetch "$acc" -O /DATABIG/sara.sebestova/SRAs/B_atticus/; done < atticus_acc.tsv
```
Then we used **fasterq-dump** in a loop to extract fastq files for every folder separately. All fastq files are in separated folders `fastq_at, fastq_gr, fastq_ros`.
Example of used command: 
```
for acc in ./*/*.sra; do fasterq-dump $acc; done
```
Before running dnaPipeTE, we need to build a library (specialized BLAST database from the input genome FASTA file). There is currently only avaliable a genome for *Bacillus rossius*.
Used command:
```
BuildDatabase -name Bacillus_r_db GCF_032445375.1_Brsri_v3_genomic.fna && RepeatModeler -database Bacillus_r_db -threads 20 -LTRStruct > Bacillus_repeatmodeler.log
```
All `SRR` folders were put in `SRR_folders` for every species, so it is less chaotic.

Before actually running dnaPipeTE, we also need to exclude mtDNA reads, because it could possibly provide fake positives if not cleaned. In `B_rossius` and `B_atticus` folder, we downloaded mitochondrion, partial genome (as for B. grandii, there is none on NCBI), named `B_ros_mt.fna` and `B_at_mt.fna`.
