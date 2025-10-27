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
Before running dnaPipeTE, we need to build a library (specialized BLAST database from the input genome FASTA file). There is currently only avaliable genome for *Bacillus rossius*.
Used command:
```
BuildDatabase -name Bacillus_r_db GCF_032445375.1_Brsri_v3_genomic.fna && RepeatModeler -database Bacillus_r_db -threads 20 -LTRStruct > Bacillus_repeatmodeler.log
```
All `SRR` folders were put in `SRR_folders` for every species, so it is less chaotic.

Before actually running dnaPipeTE, we also needed to **exclude mtDNA reads**, because it could possibly provide fake positives if not cleaned. In species folders, we downloaded mitochondrion, partial genome, named `B_ros_mt.fna`, `B_at_mt.fna`, `B_gr_mt.fa`.

We needed to index downloaded mtDNA genomes using `bwa index` and align the reads using `bwa mem` resulting in `sam` files with alignments (unmapped and mapped reads). Next it was needed to convert these files into `bam` files with only non-mitochondrial reads using `samtools` and convert them into FASTQ files `.clean.R1.fastq, .clean.R2.fastq` already cleaned files using `bedtools bamtofastq`. Example of used commands in a loop:
```
DATA_DIR="/DATABIG/sara.sebestova/SRAs/B_atticus/fastq_at"
OUT_DIR="/DATABIG/sara.sebestova/SRAs/B_atticus/cleaned_fastq"

for R1 in "$DATA_DIR"/*_1.fastq; do
    SAMPLE=$(basename "$R1" _1.fastq)
    R2="$DATA_DIR/${SAMPLE}_2.fastq"
    bwa mem B_at_mt.fna "$R1" "$R2" > "${OUT_DIR}/${SAMPLE}.mt.sam"
    samtools view -b -f 4 "${OUT_DIR}/${SAMPLE}.mt.sam" | \
    samtools sort -n -o "${OUT_DIR}/${SAMPLE}.mt.sorted.bam"
    bedtools bamtofastq -i "${OUT_DIR}/${SAMPLE}.mt.sorted.bam" \
      -fq "${OUT_DIR}/${SAMPLE}.clean.R1.fastq" \
      -fq2 "${OUT_DIR}/${SAMPLE}.clean.R2.fastq"
done
```
We blasted consensus sequences from the repeat library of *B. rossius* against mitochondrial genome of *B. rossius* to detect if any consensus sequences have significant similarity to mtDNA, so we can later filter them out of the repeat library. `consensi.fa` file is located in `/DATABIG/sara.sebestova/SRAs/genome_db/RM_3723322.TueSep301514292025$` and for blasting we used these commands: 
```
makeblastdb -in B_ros_mt.fna -dbtype nucl -out mt_db
blastn -db mt_db -query consensi.fa -evalue 0.05 -outfmt 7 > consensi_vs_mt_blast_05.txt
```
Results of this BLAST is saved in `consensi_vs_mt_blast_05.txt`. Files `mt_db` are database from `B_ros_mt.fna`. Two consensus sequences from the repeat library have strong hits to the *B. rossius* mt genome (GU001956.1). The first hit is for `rnd-5_family-55` with 95.6 % identity over 1882 bp, the second hit is for `rnd-5_family-1658` with 91.7 % identity over 1295 bp. 

In a folder `copy_consensi` we copied all `consensi.fa` files before filtering out consensus sequences with significant hits from the repeat library, just in case we made a mistake and needed to do it again. 

To know, what to use for filtering out consensus sequences, we saved relevant information about these hits present in all `consensi.fa` files to `rnd5_1658_family_list.txt` and `rnd5_1658_family_list.txt` files using grep commands: 
```
grep "^>rnd-5_family-1658" consensi* > rnd5_1658_family_list.txt
grep -E "^>rnd-5_family-55([^0-9]|$)" consensi.* > rnd5_55_family_list.txt
```
These hits were later filtered out of all `consensi.fa` files using `awk` and saved under names `filtered_consensi_*`:
```
awk '
  /^>rnd-5_family-1658/ || /^>rnd-5_family-55([^0-9]|$)/ {skip=1}
  /^>/ && !/^>rnd-5_family-1658/ && !/^>rnd-5_family-55([^0-9]|$)/ {skip=0}
  !skip
' consensi.fa > filtered_consensi_1658_55.fa

awk '
  /^>rnd-5_family-55#rRNA/ || /^>rnd-5_family-1658#Unknown/ {skip=1}
  /^>/ && !/^>rnd-5_family-55#rRNA/ && !/^>rnd-5_family-1658#Unknown/ {skip=0}
  !skip
' consensi.fa.classified > filtered_consensi_55_1658.fa.classified
```

We are running **dnaPipeTE** on every fastq file of *B. rossius* separately through the Apptainer, using only single-end R1 clean fastq files saved in `/DATABIG/sara.sebestova/SRAs/B_rossius/cleaned_fastq`. To get into the container, we have to first go into folder, where we saved `dnapipete.sif` through which we can acces the Apptainer. We always mount our directories first (directory with our clean raw reads and directory, where the repeat library is located), so we can acces it later in the container. 

```
apptainer shell \
 --bind /DATABIG/sara.sebestova/SRAs/B_rossius/cleaned_fastq:/data/fastq_parent \
 --bind /DATABIG/sara.sebestova/SRAs/genome_db/RM_3723322.TueSep301514292025:/data/repeat_lib \
 dnapipete.sif
```

Once we mounted it and we got into the container, we can run **dnaPipeTE** inside the container using command with pathway to our repeat library of *B. rossius*:
````
python3 /opt/dnaPipeTE/dnaPipeTE.py \
  -input /data/fastq_parent/SRR10323852.clean.R1.fastq \
  -output /data/fastq_parent/01_dnapipete_out_SRR10323852 \
  -cpu 2 \
  -genome_size 2000000000 \
  -genome_coverage 0.1 \
  -species B_rossius \
  -contig_length 150 \
  -RM_lib /data/repeat_lib/filtered_consensi_55_1658.fa.classified
````
Our results of **dnaPipeTE** are located in `/DATABIG/sara.sebestova/SRAs/B_rossius/cleaned_fastq` in additional folders `SRR103238*_01_dnapipete_out` where * corresponds to the specific SRR number and 01 means used genome coverage (later on, we will also run it again increasing the genome coverage). Once we have results in our `dnapipete_out` folders, we can use **dnaPT_utils** inside the Apptainer to get additional charts.
````
apptainer shell --bind /DATABIG/sara.sebestova/SRAs/B_rossius/cleaned_fastq/01_dnapipete_out_SRR10323852:/da
ta/output:rw dnapipete.sif
````
Where `:rw` means, that the folder is writable,
````
cd /data/output/

bash /opt/dnaPT_utils/dnaPT_charts.sh -I /data/output
bash /opt/dnaPT_utils/dnaPT_landscapes.sh -I /data/output
````
`dnaPT_charts` created `dnaPipeTE_charts.pdf` and `RPlots.pdf`, `dnaPT_landscapes.sh` created `dnaPipeTE_landscapes_subclass.pdf`, all in `dnapipete_out` folders.
