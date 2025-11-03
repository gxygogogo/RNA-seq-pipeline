#!/bin/bash

for PREFIX in ctrl-1 ctrl-2 ctrl-3 l-1 l-2 l-3 n-1 n-2 n-3 t-1 t-2 t-3
do
mkdir /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
cd /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
echo "bash /public1/xinyu/RNA-seq/scr/run.RNAseq.hisat2.sh ${PREFIX} mm10 10" | qsub -d ./ -N ${PREFIX} -l nodes=node03:ppn=3
done

