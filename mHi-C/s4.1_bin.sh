#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu
## Update: May 2018

################################################################
# step 4 remove duplicates among Valid Interactions and binning#
################################################################
validP=$1
validI=$2
bin=$3
dir=${validI%/*}

splitByChrom=${4:-1}
saveSplitContact=${5:-0};shift 15
chrList=("$@")


## **********************
## 4.1 remove duplicates
## **********************

if [ ! -d $dir/sorttmp ]; then
    mkdir -p $dir/sorttmp
fi

number='^[0-9]+$'
if [ "$splitByChrom" -eq 1 ];then
    ## split validPairs
    for c in ${chrList[@]} ##$(seq 1 22) X Y
    do
	if [[ $c =~ $number ]]; then
	    chrom="chr"$c
	else
	    chrom=$c
	fi

	for type in UNI MULTI
	do
	    ## remove PCR duplicates based on alignment chrom + position
	    awk -v chrom="$chrom" -v type="$type" '$2 == chrom && $13 == type {print $0}' $validP |  sort -k2,2V -k3,3n -k4 -k7,7V -k8,8n -k9 -k1 -T $dir/sorttmp | awk -v OFS="\t" 'BEGIN{c1=0;c2=0;p1=0;p2=0;s1=0;s2=0;}(c1!=$2 || c2!=$7 || p1!=$3 || p2!=$8 || s1!=$4 || s2!=$9 ){print;c1=$2;c2=$7;p1=$3;p2=$8;s1=$4;s2=$9}' >$validI.$chrom.$type
	done

	## check if Multi have duplicate alignment with Uni
	awk 'FNR==NR{a[$2$3$4$7$8$9];next};($2$3$4$7$8$9 in a)' $validI.$chrom.UNI $validI.$chrom.MULTI | awk '{print $1}' >$validI.$chrom.MultiRMdupList

    done

    ## merge all the multi-reads that overlap with uni-reads
    cat $validI.*.MultiRMdupList | sort | uniq >$validI.MultiRMdupList
    cat $validI.*.MULTI >$validI.MULTI
    ## get non-duplicated uni-reads binPairs
    cat $validI.*.UNI >$validI.UNI

    if [ "$saveSplitContact" -eq "0" ]; then
	rm -rf $validI.*.UNI
	rm -rf $validI.*.MULTI
	rm -rf $validI.*.MultiRMdupList
    fi

else
    for type in UNI MULTI
    do
	awk -v type="$type" '$13 == type {print $0}' $validP |  sort -k2,2V -k3,3n -k4 -k7,7V -k8,8n -k9 -k1 -T $dir/sorttmp | awk -v OFS="\t" 'BEGIN{c1=0;c2=0;p1=0;p2=0;s1=0;s2=0;}(c1!=$2 || c2!=$7 || p1!=$3 || p2!=$8 || s1!=$4 || s2!=$9 ){print;c1=$2;c2=$7;p1=$3;p2=$8;s1=$4;s2=$9}' >$validI.$type
    done

    ## check if Multi have duplicate alignment with Uni
    awk 'FNR==NR{a[$2$3$4$7$8$9];next};($2$3$4$7$8$9 in a)' $validI.UNI $validI.MULTI | awk '{print $1}' >$validI.MultiRMdupList

fi

## remove all the multi-reads with all their potential alignment positions if they overlap with uni-reads
awk 'FNR==NR{a[$1];next};!($1 in a)' $validI.MultiRMdupList $validI.MULTI >$validI.MULTI.noOverlapUni
rm -rf $validI.MULTI
mv $validI.MULTI.noOverlapUni $validI.MULTI

awk -v OFS="\t" '{print $1, $2, $6, $7, $11}' $validI.UNI >$validI.UNI.binPair

rm -rf $dir/sorttmp
