#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH --partition=slim18

#######################################################
## Script to call each step of mHi-C
## From the mHiC github (Ye Zheng), update May 2018
## Modified for T. brucei by Benedikt G Brink using Snakemake
## Adapted from it by Claudia Rabuffo, LMU Munich
## in collaboration with Markus Schmidt
#######################################################
######## requires: input folder with folder containing "chromosome_list.txt" (space separated contigs names) as well as "raw_reads > SAMPLE_NAME > .fq.gz files"

main(){
	set_variables
	create_sample_names
	create_folders

	conda activate mhic_env
	generate_genome_file_without_unitigs
	conda deactivate

	conda activate python2
	generate_genome_size_file
	generate_chromosome_list
	conda deactivate

	conda activate mhic_env
	faidx
	conda deactivate

	conda activate bedtools	
	generate_digestion_file_microC
	conda deactivate

	conda activate python2
	compile_cutsite
	conda deactivate

	conda activate mhic_env
	bwa_indexing
	alignment
	read_ends_pairing
	valid_fragment_filtering
	duplicates_removal
	conda deactivate

	conda activate mhic_env
	norm_bin
	generative_model_prior
	assign_multi_reads
	post_processing

	conda deactivate
	conda activate cool_env
	ploidy_correction
	merge_contigs
	cool_conversion
	cool_merge
	conda deactivate


}



###############################################################################
# Setting of global variables
###############################################################################

set_variables(){
	##############################   General   ##############
	PROJECT_PATH=$(pwd)
	ENV_PATH=/path/to/conda/environments/
    INPUT_FOLDER=input
    OUTPUT_FOLDER=output

    GENOME_VERSION=HGAP3_Tb427_2024
    GENOME_SOURCE=/path/to/genome/fasta/
    GENOME_FILE=Tb427v12

    ##############################   Preprocessing   ##############
 	ENZYME="MNase"		##MICRO C 
 	CUTSITE=0    		##MICRO C

    ##############################   s1/2   #########################
   	BWA_DIR=/path/to/conda/environments/envs/mhic_env/bin
    SAMTOOLS_DIR=/path/to/conda/environments/envs/mhic_env/bin
 	bin=${PROJECT_PATH}/bin
 	RESOLUTION=(5000)
 	SEQ_LENGTH=25
 	SAVE_FILES=0

    ##############################   s3   #########################
    REFRAG=${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_${ENZYME}.bed 		#MicroC
    refragL=0
	refragU=500

    ##############################   s4.1   #########################
	splitByChrom=0 
	saveSplitContact=0
	THREADS=64
	KMER=76

	PREFIX=$(echo ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs_mappability)

    ##############################   s4.2   #########################
	MIN_COUNT=1 #min contact counts allowed
	normMethod="None" #1. "ICE" 2. "KR" 3."None"
	# ICEminMap=0.25 ## min mappability threshold for ICE method
	# ICEmaxIter=100 ## maximum number of iteration for ICE method
	# KRsparsePerc=5 ## remove *% of sparse regions for KR method


    ##############################   s5   #########################
    splineBin=200
	priorName="uniPrior"

    ##############################   s6   #########################
	THRESHOLD=0.6

    ##############################   Programs    ##################
    BEDTOOLS_EXE=bedtools
    SAMTOOLS_EXE=samtools
    HICSD_EXE=hicsd
    NCORES=8
}

###############################################################################
# Prelude
###############################################################################

create_sample_names(){
	read -a sample_list <<< $(ls ${INPUT_FOLDER}/raw_reads)
	# echo "${sample_list[@]}"
}

create_replicate_names_for_sample(){
	read -a replicate_list <<< $(ls ${INPUT_FOLDER}/raw_reads/$1/*_R1.fq.gz | sed -e "s/_R1.fq.gz//" | sed 's/.*_//')
	# echo "${replicate_list[@]}"
}

create_folders(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
	create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				mkdir -p \
				${OUTPUT_FOLDER} \
				${OUTPUT_FOLDER}/Tbrucei_files \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s1 \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s2 \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s3_${RES} \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES} \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s5_${RES} \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES} \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}
			done
		done
	done
}



###############################################################################
# Preprocessing
###############################################################################

generate_genome_file_without_unitigs(){
	echo test && \
	seqkit seq -w 0 ${GENOME_SOURCE}/${GENOME_FILE}.fasta | grep -A 1 -E "^>BES|^>Chr" \
	> ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs.fa
 }

generate_genome_size_file(){
    python2 bin/generate_genome_size_file.py \
	   ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs.fa \
    	   > ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}.sizes
}

generate_chromosome_list(){
	read -d " " \
	-a chromosome_list < ${INPUT_FOLDER}/chromosome_list.txt
}

faidx(){
    samtools faidx ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs.fa
}

generate_digestion_file_microC(){
	${BEDTOOLS_EXE} makewindows \
	-g ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}.sizes \
	-w 100 | \
	awk 'BEGIN{ FS = OFS = "\t" } {print $1,$2,$3,$1_$2}' \
	> ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_${ENZYME}.bed
}


###############################################################################
# step 1: Trim chimeric reads, index and alignment
###############################################################################

compile_cutsite(){
	echo "Start step 1 - compile cutsite!"
	g++ -std=c++0x -o $bin/cutsite_trimming_mHiC $bin/cutsite_trimming_mHiC.cpp
}

bwa_indexing(){
	echo "Start step 1 - indexing!"
	bwa index ${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs.fa
}


alignment(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			echo "Start step 1 - alignment!"
			bash s1_bwaAlignment.sh \
				"${SAMPLE_NAME}_${REPLICATE}" \
				"${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}_without_unitigs.fa" \
				"${BWA_DIR}" \
				"${SAMTOOLS_DIR}"\
				"${INPUT_FOLDER}/raw_reads/${SAMPLE_NAME}" \
				"${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s1" \
				"${bin}" \
				"${NCORES}" \
				"${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s1/mHiC.summary_s1" \
				"${SAVE_FILES}" \
				"${SEQ_LENGTH}" \
				"${CUTSITE[@]}" &
		done
	done
	wait
}


## **************************
## step 2: Read ends pairing
## **************************


read_ends_pairing(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			echo "Start step 2 - joining read ends!"
			/path/to/conda/environments/envs/mhic_env/bin/python s2_joinEnd.py \
				-r1 $OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s1/${SAMPLE_NAME}_${REPLICATE}_1.bam \
				-r2 $OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s1/${SAMPLE_NAME}_${REPLICATE}_2.bam \
				-o ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s2/${SAMPLE_NAME}_${REPLICATE}.bam \
				-sf ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s2/mHiC.summary_s2 &
		done
	done
	wait
}


## *********************************
## step 3: Valid fragment filtering
## *********************************

valid_fragment_filtering(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Start step 3 - categorize read pairs!"
				python3 s3_categorizePairs.py \
					-f ${REFRAG} \
					-r ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s2/${SAMPLE_NAME}_${REPLICATE}.bam \
					-o ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s3_${RES} \
					-l $refragL \
					-u $refragU \
					-d ${RES} \
					-m "window" \
					-b ${RES} \
					-sf ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s3_${RES}/mHiC.summary_w${RES}_s3 &
			done
		done
	done
	wait
}



## ***************************************
## step 4.1 - Remove duplicates and binning.
## ***************************************


duplicates_removal(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Start step 4.1 - duplicates removal and binning!"
				bash s4.1_bin.sh \
						${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s3_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs \
			            ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs \
			            "$bin" \
			            "$splitByChrom" \
			            "$saveSplitContact" \
			            "${chromosome_list[@]}" &
			done
		done
	done
	wait
}


## ***************************************
## step 4.2 - Normalization and binning
## ***************************************

norm_bin(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Start step 4.2 - normalization!"
			    bash s4.2_normalization.sh \
			    	"ignored" \
			        "${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs" \
			        "$bin" \
			        "${RES}" \
			        "${MIN_COUNT}" \
			        "$normMethod" \
			        "whole" \
			        "${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}.sizes" \
			        "$KRsparsePerc" \
			        "${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/mHiC.summary_w${RES}_s4" \
			        "$splitByChrom" \
			        "$saveSplitContact" \
			        "${chromosome_list[@]}" &
			done
		done
	done
	wait
 }

## **********************
## step 5 - Build prior.
## **********************

generative_model_prior(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Starts step 5 - prior construction based on uni-reads only!"
				python3 $bin/createFitHiCFragments-fixedsize.py \
				--chrLens "${OUTPUT_FOLDER}/Tbrucei_files/${GENOME_FILE}.sizes" \
				--resolution "${RES}" \
				--outFile "$OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s5_${RES}/${SAMPLE_NAME}_${REPLICATE}_${RES}.uni.fragments.mHiC"

				python3 s5_prior.py \
				-f $OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s5_${RES}/${SAMPLE_NAME}_${REPLICATE}_${RES}.uni.fragments.mHiC \
				-i ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPairCount.uni \
				-o $OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s5_${RES} \
				-b $splineBin \
				-L $((${RES} * 2)) \
				-r ${RES} \
				-p 1 &
			done
		done
	done
	wait
}



## ************************************************************************************
## step 6.1 - Generative model to assign probability to multi-reads potential alignments.
## ************************************************************************************


assign_multi_reads(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Starts step 6.1 - assign probability to multi reads alignments"
				awk -v OFS="=" '{print $2, $3, $4, $5}' \
					${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.MULTI.binPair.multi \
					| sort -u \
					>${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.MULTI.binPair.multiKeys

				python s6_em.py \
					-p $OUTPUT_FOLDER/${SAMPLE_NAME}_${REPLICATE}/s5_${RES}/s5_prior.mhic \
					-u ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPairCount.uni \
					-m ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.MULTI.binPair.multi  \
					-mk ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.MULTI.binPair.multiKeys \
					-t ${THRESHOLD} \
					-o "${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}" \
					-f ${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.multi &
			done
		done
	done
	wait
}



## ************************************************************************************
## step 6.2 - Post-mHiC processing.
## ************************************************************************************

post_processing(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				echo "Starts step 6.2 - post mHiC processing"
				awk -v OFS="\t" \
				-v fT=$THRESHOLD '$6>fT {print $2, $3, $4, $5}' \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.multi.mHiC | \
				sort | \
				uniq -c | \
				awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' \
				>${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.multi.mHiC.binPairCount.multi ## get binPair Count for multi-reads
					
				#cat $uni $multiOut.binPairCount.multi
				cat ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s4_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPairCount.uni \
				${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.multi.mHiC.binPairCount.multi \
				| sort -k1,1V -k2,2n -k3,3V -k4,4n \
				| awk -v OFS="\t" '{a[$1" "$2" "$3" "$4]+=$5}END{for (i in a) print i,a[i]}' \
				| sort -k1,1V -k2,2n -k3,3V -k4,4n \
				>${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.uniMulti ## merged with uni-reads binpair count
			done
		done
	done
}


############ PLOIDY CORRECTION AND CONVERSION TO ENTIRE CHROMOSOMES
#### THIS STEP (parameter s and name of output) CHANGES WITH RESOLUTION - MANUALLY CORRECT IT! s should be half the resolution
ploidy_correction(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				python3 ${PROJECT_PATH}/bin/ploidy_normalization_matrix.py \
						${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.validPairs.binPair.uniMulti \
						> ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.guided_ploidy.uniMulti 
				cat ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_5000_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.guided_ploidy.uniMulti  | \
						awk -v s=2500 '{print $1,$2-s,$2+s,$3,$4-s, $4+s,$5}' | tr " " "\t" | gzip > ${OUTPUT_FOLDER}/cool_files/${SAMPLE_NAME}_${REPLICATE}_5000.guided_ploidy.bg2.gz &
			done
		done
	done
	wait
}

#### THIS STEP CHANGES WITH RESOLUTION - MANUALLY CORRECT IT! s should be half the resolution
merge_contigs(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			for RES in ${RESOLUTION[@]} 
			do
				python3 ${PROJECT_PATH}/bin/merge_contigs_in_uniMulti_file.py \
						/path/to/ctsizes/file.ctsizes \
						${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.guided_ploidy.uniMulti \
						>  ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_${RES}_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.entire_chrs.uniMulti
				cat ${OUTPUT_FOLDER}/${SAMPLE_NAME}_${REPLICATE}/s6_5000_${THRESHOLD}/${SAMPLE_NAME}_${REPLICATE}.entire_chrs.uniMulti  | \
						awk -v s=2500 '{print $1,$2-s,$2+s,$3,$4-s, $4+s,$5}' | tr " " "\t" | gzip > ${OUTPUT_FOLDER}/cool_files/${SAMPLE_NAME}_${REPLICATE}_5000.entire_chrs.bg2.gz &
			done
		done
	done
	wait
}


cool_conversion(){
	for SAMPLE_NAME in ${sample_list[@]}
	do
		create_replicate_names_for_sample ${SAMPLE_NAME}
		for REPLICATE in ${replicate_list[@]}
		do
			#### entire chromosomes
			cooler load -f bg2 /path/to/entire_chromosome/files/file_entire_chromosome.sizes:5000 \
					${OUTPUT_FOLDER}/cool_files/${SAMPLE_NAME}_${REPLICATE}_5000.entire_chrs.bg2.gz \
					${OUTPUT_FOLDER}/cool_files/mhic_${SAMPLE_NAME}_${REPLICATE}_5000.entire_chrs.cool
			cooler zoomify --nproc 8 --out ${OUTPUT_FOLDER}/cool_files/mhic_${SAMPLE_NAME}_${REPLICATE}_5000.entire_chrs.mcool --resolutions 500000,100000,50000,25000,10000,5000 --balance \
					${OUTPUT_FOLDER}/cool_files/mhic_${SAMPLE_NAME}_${REPLICATE}_5000.entire_chrs.cool
		done
	done
	wait
}

#### sample names of samples that have to be merged should be put manually 
cool_merge(){
	#### entire chromosomes
	cooler merge ${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_merged_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib1_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib2_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib3_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib4_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib5_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib6_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib7_5000.entire_chrs.cool \
					${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_Lib8_5000.entire_chrs.cool
	cooler zoomify --nproc 8 --out ${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_merged_5000.entire_chrs.mcool --resolutions 500000,100000,50000,25000,10000,5000 --balance ${OUTPUT_FOLDER}/cool_files/{SAMPLE_NAME}_merged_5000.entire_chrs.cool	
}


main