#!/bin/bash

# -q 19: 质量值低于此值认为是低质量碱基，默认15，表示质量为Q15
# -u 50: 限制低质量碱基的百分比，默认40，表示40%
# -n 5: 限制N的数量
# --length_required 150: 最小长度值
# --length_limit 150: 最大长度值
alias fastp=/home/xinyu/Software/fastp

for PREFIX in ctrl-1 ctrl-2 ctrl-3 l-1 l-2 l-3 n-1 n-2 n-3 t-1 t-2 t-3
do
echo ${PREFIX}
mkdir /public1/xinyu/RNA-seq/CleanData/${PREFIX}
cd /public1/xinyu/RNA-seq/CleanData/${PREFIX}
fastp \
   -i /public1/xinyu/RNA-seq/RawData/${PREFIX}/${PREFIX}_1.fq.gz \
   -o /public1/xinyu/RNA-seq/CleanData/${PREFIX}/${PREFIX}_1.fq.gz \
   -I /public1/xinyu/RNA-seq/RawData/${PREFIX}/${PREFIX}_2.fq.gz \
   -O /public1/xinyu/RNA-seq/CleanData/${PREFIX}/${PREFIX}_2.fq.gz \
   -q 19 \
   -u 50 \
   -n 5 \
   --length_required 150 \
   --length_limit 150
mv *.html ${PREFIX}.html
mv *.json ${PREFIX}.json
done

