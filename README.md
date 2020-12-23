# ATAC-seq

```
##1.fastq QC
fastqc read1 read2

###2.mapping

bwa mem -t 20 -M ref.fa read1.fastq read2.fastq > sample.sam

samtools view -bS sample.sam -o sample.bam
 
samtools flagstat sample.bam

samtools view -f 3 -q 30 -o sample_q30.bam sample.bam # -f 3 specifies only properly-paired reads

samtools sort -@ 20 -o sample_q30_sort.bam sample_q30.bam

java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=sample_q30_sort.bam O=sample_q30_sort_dedup.bam M=sample_sort.markduplicates.log

samtools index sample_q30_sort_dedup.bam


ataqv --name sample --autosomal-reference-file "file name" --metrics-file sample.ataqv.json.gz rice sample_q30_sort_dedup.bam > sample.ataqv.out
mkarv my_fantastic_experiment sample.ataqv.json.gz

# in R

#### fragSizeDist
library(ATACseqQC)
bamfile <- sample_q30_sort_dedup.bam

fragSizeDist(bamfiles, bamfiles.labels=bbamfile.sample.ID, index = bamfiles, ylim = NULL, logYlim = NULL )

pdf(paste0(bamfile.sample.ID, ".fragment.size.distribution.pdf"), width =10, height=8) 
fragSize <- fragSizeDist(bamFiles=bamfile, bamFiles.labels=bamfile.sample.ID)
dev.off()


# Estimate library complexity
# Define sample files
treat1 <- "treat1.sorted.filt.noMT.bam"
treat1.bai <- "treat1.sorted.filt.noMT.bam.bai"

# Calculate duplication frequency matrix
treat1.dups <- readsDupFreq(treat1, index=treat1.bai)

# Estimate library complexity
treat1.complexity <- estimateLibComplexity(treat1.dups, times=100, interpolate.sample.sizes=seq(0.1, 1, by=0.01)


##
0x2	PROPER_PAIR	 每对短序列都被 aligner 合适的比对上了


####Tn5 reads shift

###bam文件转成bw文件 IGV可视化
bamCoverage --bam a.bam -o a.SeqDepthNorm.bw --binSize 1 --normalizeUsing RPKM --outFileFormat bedgraph 

## TSS/TES enrichment
computeMatrix reference-point \ # choose the mode
       --referencePoint TSS \ # alternatives: TSS, TES, center
       -b 3000 -a 3000 \ # define the region you are interested in
       -R genes.gtf \
       -S testFile.bw  \
       --skipZeros \
       -o matrix1_TSS.gz \ # to be used with plotHeatmap and plotProfile
       --outFileSortedRegions regions1_genes.bed 

      plotHeatmap -m matrix.mat.gz \
      -out ExampleHeatmap1.png \



#call peaks

#rice 的有效基因组大小3.0e8
macs2 callpeak -t *.bam -f bam --shift 75 --extsize 150 --nomodel --keep-dup all -B --SPMR -g 3.0e8 --outdir Macs2_out

bam2bed

macs2 callpeak -t ATAC1.rmdup.bed -n ATAC1 -g 1.0e8 --call-summit -f BAMPE -nomodel -B --SPMR --extsize 200 

macs2 callpeak -t treat1_sub.minimal.bedpe -f BEDPE -n ATAC1 -g 3.0e8 --broad --broad-cutoff 0.05 --keep-dup all


# Call significant broad peaks (FDR < 0.05) on each individual replicate
macs2 callpeak -t treat1_sub.minimal.bedpe -f BEDPE -n treat1 --broad --broad-cutoff 0.05 --keep-dup all

macs2 callpeak --nomodel -t *.bam -q 0.01 -f BAM  --keep-dup all




######  HMMRATAC  Quick Start
samtools sort ATACseq.bam -o ATACseq.sorted.bam
samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai
samtools view -H ATACseq.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info
java -jar HMMRATAC_V1.2.4_exe.jar -b ATACseq.sorted.bam -i ATACseq.sorted.bam.bai -g genome.info



Filter HMMRATAC output by the score, if desired. Score threshold will depend on dataset, score type and user preference. A threshold of 10 would be:

awk -v OFS="\t" '$13>=10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak

To filter the summit file by the same threshold:

awk -v OFS="\t" '$5>=10 {print}' NAME_summits.bed > NAME.filteredSummits.bed

## NOTE: HMMRATAC will report all peaks that match the structure defined by the model, including weak peaks. 
## Filtering by score can be used to retain stronger peaks. Lower score = higher sensitivity and lower precision, Higher score = lower sensitivity and higher precision.


### 在输入HMMRATAC之前，可以删除重复或低Mapping质量的读数，尽管HMMRATAC默认会执行这些功能。默认情况下，HMMRATAC会删除重复的读数（而且这个功能目前是硬编码的，无法关闭）。
### HMMRATAC还将删除那些MapQ（映射质量分数）低于30的读数。这个值可以通过命令行选项（使用-q或--minmapq选项）改变。



## Peak annotation
homer

wget http://homer.salk.edu/homer/configureHomer.pl  #首先下载configureHomer.pl
perl configureHomer.pl -install  #使用configureHomer.pl配置Homer

perl configureHomer.pl  -install rice.IRGSP-1.0  ##准备参考基因组的注释信息

annotatePeaks.pl peak.bed rice.IRGSP-1.0 > peak.annotation.xls
```
