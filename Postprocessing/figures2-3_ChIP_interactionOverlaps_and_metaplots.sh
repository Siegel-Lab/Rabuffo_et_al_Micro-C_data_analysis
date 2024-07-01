#!/bin/bash

############################################
###### 	loop calling with Mustache	   #####
############################################
conda activate mustache

./mustache/mustache/mustache.py  -f /path/to/file.mcool -r 5kb -o output_called_loops_5kb.tsv -pt 0.02 -st 0.85

conda deactivate 

############################################
###### 		ChIP seq data analysis 	   #####
############################################

module load ngs/samtools/1.9
module load ngs/MACS2/2.1.2
module load ngs/bedtools2/2.28.0
module load ngs/UCSCutils/3.4.1
module load ngs/bwa/0.7.16 
module load ngs/samtools/1.8
module load ngs/deeptools/3.5.0
module load gcc/9.2.0

bwa index /path/to/genome/Tb427v12.fasta

##### H2AZ
bwa mem -t 16 /path/to/genome/Tb427v12.fasta H2AZ_ChIP_file.fastq.gz | samtools view -bh > H2AZ_ChIP_output.bam
samtools sort -@ 16 H2AZ_ChIP_output.bam -o H2AZ_ChIP_output.sorted.bam
samtools markdup -r H2AZ_ChIP_output.sorted.bam H2AZ_ChIP_output_nodups.sorted.bam
samtools index H2AZ_ChIP_output_nodups.sorted.bam
bamCoverage -b H2AZ_ChIP_output_nodups.sorted.bam -o H2AZ_ChIP_output_nodups.bw

macs2 callpeak -t ./bamfiles/SRR5466320_entire_chrs_nodups.sorted.bam -f BAM -g 3.5e+7 --max-gap 500 -n H2AZ_nm_extsize300_nodups_q001_maxgap500 -q 0.01 --broad --nomodel --extsize 300 --cutoff-analysis --trackline --outdir ./macs2

##### SCC1
bwa mem -t 16 /path/to/genome/Tb427v12.fasta SCC1_ChIP_file.fastq.gz | samtools view -bh > SCC1_ChIP_output.bam
samtools sort -@ 16 SCC1_ChIP_output.bam -o SCC1_ChIP_output.sorted.bam
samtools index SCC1_ChIP_output.sorted.bam

bwa mem -t 16 /path/to/genome/Tb427v12.fasta SCC1_input_file.fastq.gz | samtools view -bh > SCC1_input_output.bam
samtools sort -@ 16 SCC1_input_output.bam -o SCC1_input_output.sorted.bam
samtools index SCC1_input_output.sorted.bam

bamCompare -b1 SCC1_ChIP_output.sorted.bam \
 			 -b2 SCC1_input_output.sorted.bam \
 			 --outFileFormat bigwig \
 			 -o SCC1_ratio.bw \
 			 --operation ratio \
 			 --ignoreDuplicates

macs2 callpeak -t SCC1_ChIP_output.sorted.bam -c SCC1_input_output.sorted.bam -f BAM -g 3.5e+7 -n SCC1 --nomodel --extsize 200 --trackline --outdir ./macs2


#############################################################
###### overlapping loops with annotations or peak files #####
#############################################################
###### 
# # usage: python3 overlapping_loop_annotation_counts.py <loop_input_file> <require_peaks_on_both_axes> (note: true / false) <annotation_input_file> <annotation_type>
#### note: the annotation file should have the third column with the type of feature, in the following example TSS
#### this is true also for called peaks: the gff file should have the third column with the type of annotation (e.g., peak)
# ## both loops
python3 bin/overlapping_loop_annotation_counts.py \
		output_called_loops_5kb.tsv  \
		TRUE \
 		/path/to/called_peaks.gff \
 		TSS \
		> output_number_of_overlaps.txt

# ## at least one loop coord
python3 bin/overlapping_loop_annotation_counts.py \
		output_called_loops_5kb.tsv  \
		FALSE \
 		/path/to/called_peaks.gff \
 		TSS \
		>> output_number_of_overlaps.txt

############################################
###### 			Metaplots			   #####
############################################

### plot along the annotated features
computeMatrix scale-regions -S /path/to/ChIP_over_input.bw  -R /path/to/annotation_file.bed -a 0 -b 0 -o output_matrix.mat.gz
plotHeatmap -m output_matrix.mat.gz --heatmapWidth 5 -out output_plot.pdf

#### plot centered around the annotated features
computeMatrix scale-regions -S /path/to/ChIP_over_input.bw  -R /path/to/annotation_file.bed --referencePoint center -a 15000 -b 15000  -o output_matrix.mat.gz
plotHeatmap -m output_matrix.mat.gz --heatmapWidth 5 -out output_plot.pdf
