#!/bin/bash

WORKDIR="/path/to/you/project"
RAW_READS_DIR="$WORKDIR/raw_reads"
CLEANED_DIR="$WORKDIR/cleaned_reads"
QC_DIR="$WORKDIR/qc_reports"
THREADS=20
ASSEMBLY_DIR="$WORKDIR/idba_assembly"
MAPPING_DIR="$WORKDIR/mapping"
GENE_DIR="$WORKDIR/gene"
EGGNOG_DIR="$WORKDIR/EggNOG"
source $HOME/.bashrc
MPA_DIR="$WORKDIR/metaphlan"
HUMANN_DIR="$WORKDIR/humann"
mkdir -p $QC_DIR/initial $QC_DIR/final $CLEANED_DIR
mkdir -p $ASSEMBLY_DIR $MAPPING_DIR $ASSEMBLY_DIR/$sample
mkdir -p $GENE_DIR
mkdir -p $EGGNOG_DIR
mkdir -p $MPA_DIR
mkdir -p $HUMANN_DIR
source activate your_env_name

for sample in `cat $WORKDIR/test_list.txt`
do
prefetch $sample -O $RAW_READS_DIR
 
fastq-dump --gzip --split-3 -O $RAW_READS_DIR $RAW_READS_DIR/$sample/${sample}.sra 
 
r1=$RAW_READS_DIR/${sample}_1.fastq.gz
r2=$RAW_READS_DIR/${sample}_2.fastq.gz
c1=$CLEANED_DIR/${sample}_clean_R1.fastq.gz
c2=$CLEANED_DIR/${sample}_clean_R2.fastq.gz
fastqc -f fastq -t $THREADS ${r1} -o $QC_DIR/initial
fastqc -f fastq -t $THREADS ${r2} -o $QC_DIR/initial
bbduk.sh in=${r1} in2=${r2} out=${c1} out2=${c2} \
            qtrim=rl trimq=20 mlf=0.33 threads=$THREADS ref="/home/ug2023/ug523111910118/Biosofts/bbmap/resources/adapters.fa"
fastqc -f fastq -t $THREADS ${c1} -o $QC_DIR/final
fastqc -f fastq -t $THREADS ${c2} -o $QC_DIR/final
 
mkdir -p $ASSEMBLY_DIR/$sample
gunzip -c ${c1} > $CLEANED_DIR/${sample}_clean_R1.fastq
gunzip -c ${c2} > $CLEANED_DIR/${sample}_clean_R2.fastq
 
fq2fa --merge --filter $CLEANED_DIR/${sample}_clean_R1.fastq $CLEANED_DIR/${sample}_clean_R2.fastq $CLEANED_DIR/${sample}_clean.fa
 
idba_ud -r $CLEANED_DIR/${sample}_clean.fa -o $ASSEMBLY_DIR/$sample
 
bowtie2-build \
    --threads $THREADS \
    $ASSEMBLY_DIR/${sample}/contig.fa \
    $MAPPING_DIR/${sample}_contig_index
 
bowtie2 \
    -x $MAPPING_DIR/${sample}_contig_index \
    -1 $CLEANED_DIR/${sample}_clean_R1.fastq.gz \
    -2 $CLEANED_DIR/${sample}_clean_R2.fastq.gz \
    -S $MAPPING_DIR/${sample}_mapped.sam \
    --threads $THREADS \
    --no-unal \
    --sensitive \
    --dovetail
 
samtools view -@ $THREADS -bS $MAPPING_DIR/${sample}_mapped.sam | \
samtools sort -@ $THREADS -o $MAPPING_DIR/${sample}_mapped_sorted.bam
samtools flagstat $MAPPING_DIR/${sample}_mapped_sorted.bam > $MAPPING_DIR/${sample}_mapping_stats.txt
 
rm $CLEANED_DIR/${sample}_clean_R1.fastq
rm $CLEANED_DIR/${sample}_clean_R2.fastq
prodigal -p meta -i $ASSEMBLY_DIR/${sample}/contig.fa -a $GENE_DIR/${sample}.faa -d $GENE_DIR/${sample}.fna -o $GENE_DIR/${sample}.gff -f gff -q
 
emapper.py -i $GENE_DIR/${sample}.faa -o $EGGNOG_DIR/${sample}
 
metaphlan $CLEANED_DIR/${sample}_clean_R1.fastq.gz,$CLEANED_DIR/${sample}_clean_R2.fastq.gz \
    --input_type fastq \
    --bowtie2out $MPA_DIR/${sample}_metaphlan_bowtie2.bz2 \
    -o $MPA_DIR/${sample}_metaphlan_profile.tsv \
    --nproc $THREADS \
    --tax_lev 's'
 
cat $c1 $c2 > $CLEANED_DIR/${sample}_clean_cb.fastq.gz
humann \
    --input  $CLEANED_DIR/${sample}_clean_cb.fastq.gz \
    --output $HUMANN_DIR \
    --threads $THREADS \
    --input-format fastq.gz \
    --o-log $HUMANN_DIR/${sample}_humann3.log
 
humann_regroup_table \
    --input $HUMANN_DIR/${sample}_clean_cb_genefamilies.tsv \
    --groups uniref90_ko \
    --output $HUMANN_DIR/${sample}_ko_pathways.tsv
humann_renorm_table \
    --input $HUMANN_DIR/${sample}_ko_pathways.tsv \
    --output $HUMANN_DIR/${sample}_ko_pathways_cpm.tsv \
    --units cpm
done
