#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH --partition=slim18

##### note: fq.gz files in input folder
### Pipeline according to open2C/distiller

main(){
	set_variables
	create_folders
	module load ngs/bwa/0.7.16
	module load ngs/samtools/1.9
	
	conda activate pairtools
	mapping
	rm_tmp
	dedup
	split
	indexing
	ploidy_correction
	merge_contigs
	conda deactivate

	conda activate cool_env
	count_interactions
	cool_conversion
	merge
}


###############################################################################
# Setting of global variables
###############################################################################

set_variables(){
	##############################   General   ##############
	PROJECT_PATH=$(pwd)
	MAPQ30_FOLDER=$(pwd)/mapq30_cool

    GENOME_VERSION=Tb427v12
    GENOME_SOURCE=/path/to/genome/folder/
    GENOME_FILE=Tb427v12

    INPUT_FOLDER=input
    OUTPUT_FOLDER=output
    TEMPORARY_DIR=tmpdir

    TASK_TMP_DIR=$(mktemp -d -p ./ distiller.tmp.XXXXXXXXXX)

    NCORES=16
    NPROC=8

    RESOLUTION=1000
}

create_folders(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		READ_TWO=$(echo ${READ_FILE} | sed "s/_R1/_R2/")
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
			mkdir -p \
			${OUTPUT_FOLDER} \
			${OUTPUT_FOLDER}/${SAMPLE_NAME}  \
			${MAPQ30_FOLDER}
	done
}


mapping(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		READ_TWO=$(echo ${READ_FILE} | sed "s/_R1/_R2/")
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		touch ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.bam
		fastp -q 20 --json ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.fastp.json --html ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.fastp.html -i ${INPUT_FOLDER}/${READ_FILE} -I ${INPUT_FOLDER}/${READ_TWO} --stdout | \
		bwa mem -p -t ${NCORES} -v 3 -SP ${GENOME_SOURCE}/${GENOME_FILE}.fasta - | \
		pairtools parse --drop-sam --drop-readid --add-columns mapq,XA --walks-policy mask -c ${GENOME_SOURCE}/${GENOME_FILE}.sizes | \
		pairtools sort --nproc ${NCORES}  -o ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.pairsam.gz --tmpdir $TASK_TMP_DIR  | cat &
	done
	wait
}

rm_tmp(){
	rm -rf $TASK_TMP_DIR
}


dedup(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.pairsam.gz | \
		pairtools dedup --max-mismatch 3 --mark-dups --output \
		>( pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.nodups.pairs.gz \
			--output-sam ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.nodups.bam ) --output-unmapped \
		>( pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.unmapped.pairs.gz \
			--output-sam ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.unmapped.bam ) --output-dups \
		>( pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.dups.pairs.gz \
			--output-sam ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.dups.bam ) --output-stats ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.dedup.stats | cat &
	done
	wait
}


split(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		pairtools split ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}_dedup.pairs --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.nodups.pairs.gz &
	done
	wait
}

indexing(){
	pairix ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.nodups.pairs.gz
}

ploidy_correction(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		python3 bin/ploidy_normalization_pairs.py <(bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.nodups.pairs.gz) \
		> ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.guided_ploidy.nodups.pairs
		bgzip ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.guided_ploidy.nodups.pairs &
	done
	wait
}

merge_contigs(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		zcat ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.guided_ploidy.nodups.pairs.gz | \
		python bin/merge_contigs_in_pairs_file.py bin/Tb427v12_entire_chrs_conversion.ctsizes - | \
		gzip > ${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.entire_chrs.nodups.pairs.gz &
	done
	wait
}


cool_conversion(){
	for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
	do
		SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
		cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly ${GENOME_VERSION} /path/to/chromosome/size/file.sizes:${RESOLUTION} \
		${OUTPUT_FOLDER}/${SAMPLE_NAME}/${SAMPLE_NAME}.entire_chrs.nodups.pairs.gz \
		${MAPQ30_FOLDER}/${SAMPLE_NAME}_${RESOLUTION}.entire_chrs.cool 
	done
	wait
}

##### names of libraries to merge have to be put manually

merge(){
	cooler merge ${MAPQ30_FOLDER}/SAMPLEID_1000_entirechrs_merged.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib1_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib2_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib3_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib4_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib5_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib6_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib7_1000.entire_chrs.cool \
	${MAPQ30_FOLDER}/SAMPLEID_Lib8_1000.entire_chrs.cool
	cooler zoomify --nproc 8 --out ${MAPQ30_FOLDER}/SAMPLEID_Lib1_1000.entire_chrs.mcool --resolutions 100000,50000,25000,10000,5000,2000,1000 --balance --balance-args "--tol 0.05" ${MAPQ30_FOLDER}/SAMPLEID_Lib1_1000.entire_chrs.cool
}

main
