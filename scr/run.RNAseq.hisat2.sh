#!/bin/bash

PREFIX=${1}
GENOME_ASSEMBLY_NAME=${2}
cpu_num=${3}

output=/public2/TangLabData/ProcessedData/RNA-seq/${PREFIX}
#input=/data2/TangLabData/CleanData/RNA-seq/${PREFIX}
input=/public2/TangLabData/CleanData/RNA-seq/${PREFIX}

#fastq_1=${input}/${PREFIX}_1.fq.gz
#fastq_2=${input}/${PREFIX}_2.fq.gz

#fastq_1=${input}/${PREFIX}_R1.fastq.gz
#fastq_2=${input}/${PREFIX}_R2.fastq.gz

fastq_1=${input}/${PREFIX}_R1.fq.gz
fastq_2=${input}/${PREFIX}_R2.fq.gz

#fastq_1=${input}/${PREFIX}_R1.fq
#fastq_2=${input}/${PREFIX}_R2.fq

if [ "${GENOME_ASSEMBLY_NAME}" == "panTro5" ];
then
GENOME=/data/public/refGenome/hisat2_index/panTro5/panTro5
GENOME_SIZE=/data/public/refGenome/bwa_index/panTro5/panTro5.chrom.sizes
GENOME_ANNOTATION=/data1/tang/3D_evolution/primates_annotation/panTro5/Pan_troglodytes.Pan_tro_3.0.97.chr.gtf

elif [ "${GENOME_ASSEMBLY_NAME}" == "gorGor4" ];
then
GENOME=/data/public/refGenome/hisat2_index/gorGor4/gorGor4
GENOME_SIZE=/data/public/refGenome/bwa_index/gorGor4/gorGor4.chrom.sizes
GENOME_ANNOTATION=/data1/tang/3D_evolution/primates_annotation/gorGor4/Gorilla_gorilla.gorGor4.97.chr.gtf

elif [ "${GENOME_ASSEMBLY_NAME}" == "rheMac8" ];
then
GENOME=/data/public/refGenome/hisat2_index/rheMac8/rheMac8
GENOME_SIZE=/data/public/refGenome/bwa_index/rheMac8/rheMac8.chrom.sizes
GENOME_ANNOTATION=/data1/tang/3D_evolution/primates_annotation/rheMac8/Macaca_mulatta.Mmul_8.0.1.97.chr.gtf

elif [ "${GENOME_ASSEMBLY_NAME}" == "hg38" ];
then
GENOME=/data/public/refGenome/hisat2_index/hg38/hg38
GENOME_SIZE=/data/public/refGenome/bwa_index/hg38/ChromInfo.txt
GENOME_ANNOTATION=/data1/tang/3D_evolution/primates_annotation/hg38/Homo_sapiens.GRCh38.97_withchr.gtf

elif [ "${GENOME_ASSEMBLY_NAME}" == "mm10" ];
then
GENOME=/data/public/refGenome/hisat2_index/mm10/genome
GENOME_SIZE=/data/public/refGenome/bwa_index/mm10/mm10.chrom.sizes
GENOME_ANNOTATION=/data/public/Gene.annotation/mm10/Mus_musculus.GRCm38.98.gtf
fi

##############################################
## Tools
hisat2=/data/public/software/hisat2-2.1.0/hisat2
samtools=/data/public/software/samtools-1.3.1/samtools
bedtools=/data/public/software/bedtools.2.25.0/bin/bedtools
bedGraphToBigWig=/home/yuhan/Software/others/bin/bedGraphToBigWig

echo -e "RNA-seq \ninput: ${input}\noutput: ${output}"
PREFIX=${output}/${PREFIX}

echo "mapping"
## mapping
${hisat2} -x ${GENOME} -1 ${fastq_1} -2 ${fastq_2} -p ${cpu_num} -S ${PREFIX}.sam

## sam2bam
${samtools} view -b -o ${PREFIX}.bam  -@ ${cpu_num} - < ${PREFIX}.sam

## sort bam by coordinate
${samtools} sort --threads ${cpu_num} -o ${PREFIX}.sorted.bam ${PREFIX}.bam

## delete unsorted sam
rm -f ${PREFIX}.sam

##unique map
${samtools} view -b -o ${PREFIX}.uniq.bam -q 20 -@ ${cpu_num} ${PREFIX}.sorted.bam
${samtools} index ${PREFIX}.uniq.bam
#rm -f ${PREFIX}.sorted.bam

echo -e "counting\n"
## counting library-type
htseq-count -n ${cpu_num} -f bam -r pos -t exon --max-reads-in-buffer 1000000000 --stranded reverse -m union ${PREFIX}.uniq.bam ${GENOME_ANNOTATION} > ${PREFIX}.count 2>${PREFIX}.counting.log



## read bam convered to bed and calculate coverage
#${bedtools} bamtobed -split -i ${PREFIX}.uniq.bam | awk '$5>=20{if($1!~"[_M]") print $0;}' | LANG=C sort -k1,1 -k2,2n > ${PREFIX}.bed
#${bedtools} genomecov -bg -i ${PREFIX}.bed -g ${GENOME_SIZE} > ${PREFIX}.bedgraph
#${bedGraphToBigWig} ${PREFIX}.bedgraph ${GENOME_SIZE} ${PREFIX}.bigwig

##############################################################################################
echo -e "spliting strand-specific reads\n"

if [ -f ${PREFIX}.bed12 ]; then
    rm ${PREFIX}.bed12
fi

for index in `awk '($1!~/M/)&&($1!~/_/){print $1}' ${GENOME_SIZE}`
do
    echo "Processing chr${index}..."

    ${samtools} view -b -o ${PREFIX}.flag64_first.bam -f 64 ${PREFIX}.uniq.bam ${index}
    ${samtools} view -b -o ${PREFIX}.flag128_second.bam -f 128 ${PREFIX}.uniq.bam ${index}

    ${bedtools} bamtobed -bed12 -i ${PREFIX}.flag128_second.bam > ${PREFIX}.tmp.bed12
    ${bedtools} bamtobed -bed12 -i ${PREFIX}.flag64_first.bam | perl -lane 'if( $F[5] eq "+" ){ $F[5] = "-"; print join( "\t", @F ) } else { $F[5] = "+"; print join( "\t", @F ) }' >> ${PREFIX}.tmp.bed12

    cat ${PREFIX}.tmp.bed12 | LANG=C sort -k1,1 -k2,2n >> ${PREFIX}.bed12

    rm ${PREFIX}.tmp.bed12
    rm ${PREFIX}.flag64_first.bam
    rm ${PREFIX}.flag128_second.bam
done

echo -e "Done"
