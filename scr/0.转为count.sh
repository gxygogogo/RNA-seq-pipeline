#!/bin/bash

GENOME_ANNOTATION=/data/public/Gene.annotation/mm10/Mus_musculus.GRCm38.98.gtf

for PREFIX in ctrl-1 ctrl-2 ctrl-3 l-1 l-2 l-3 n-1 n-2 n-3 t-1 t-2 t-3
do
mkdir /public1/xinyu/RNA-seq/Count/${PREFIX}
cd /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
echo "htseq-count -f bam -r pos -t exon --max-reads-in-buffer 1000000000 --stranded reverse -m union ${PREFIX}.uniq.bam ${GENOME_ANNOTATION} > ../../Count/${PREFIX}/${PREFIX}.count 2>../../Count/${PREFIX}/${PREFIX}.counting.log" | qsub -N ${PREFIX} -l nodes=node02:ppn=1 -d ./
done

