# RNA-seq pipeline
## Step1. Raw data QC
### 使用fastp进行过滤
```{shell}
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
```
### 质控结果
![QC结果](https://github.com/gxygogogo/RNA-seq-pipeline/blob/main/img/qc.png)
## Step2. Genome mapping
### 比对到指定基因组
```{shell}
for PREFIX in ctrl-1 ctrl-2 ctrl-3 l-1 l-2 l-3 n-1 n-2 n-3 t-1 t-2 t-3
do
mkdir /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
cd /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
echo "bash /public1/xinyu/RNA-seq/scr/run.RNAseq.hisat2.sh ${PREFIX} mm10 10" | qsub -d ./ -N ${PREFIX} -l nodes=node03:ppn=3
done
```
### 转换为count
```{shell}
GENOME_ANNOTATION=/data/public/Gene.annotation/mm10/Mus_musculus.GRCm38.98.gtf

for PREFIX in ctrl-1 ctrl-2 ctrl-3 l-1 l-2 l-3 n-1 n-2 n-3 t-1 t-2 t-3
do
mkdir /public1/xinyu/RNA-seq/Count/${PREFIX}
cd /public1/xinyu/RNA-seq/ProcessData/${PREFIX}
echo "htseq-count -f bam -r pos -t exon --max-reads-in-buffer 1000000000 --stranded reverse -m union ${PREFIX}.uniq.bam ${GENOME_ANNOTATION} > ../../Count/${PREFIX}/${PREFIX}.count 2>../../Count/${PREFIX}/${PREFIX}.counting.log" | qsub -N ${PREFIX} -l nodes=node02:ppn=1 -d ./
done
```
## Step3. DE gene analysis
使用DESeq2识别差异基因并进行可视化等下游分析。<br>
![DE结果](https://github.com/gxygogogo/RNA-seq-pipeline/blob/main/img/de.png)
## Step4. Gene Enrichment
* 对差异基因进行GO、KEGG富集分析
* 对指定的通路进行GSEA分析
![GSEA结果](https://github.com/gxygogogo/RNA-seq-pipeline/blob/main/img/gsea.png)
