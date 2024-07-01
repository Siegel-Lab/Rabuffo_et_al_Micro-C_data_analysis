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
resolution=$4
minCount=${5:-1}

normMethod=${6:-"ICE"}
mappFile=$7
minMap=${8:-0.5}
maxIter=${9:-350}
normChrom=${10:-"whole"}
chromSizeFile=${11}
sparsePerc=${12:-10}
summaryFile=${13:-""}
splitByChrom=${14:-1}
saveSplitContact=${15:-0};shift 15
chrList=("$@")

if [ "$normMethod" != "ICE" ] && [ "$normMethod" != "ice" ] && [ "$normMethod" != "KR" ] && [ "$normMethod" != "kr" ]; then
    if [ "$normMethod" != "None" ] &&  [ "$normMethod" != "NONE" ] &&  [ "$normMethod" != "none" ] &&  [ "$normMethod" != "" ]; then
	echo "Normalization methods can only be either ICE or KR! Otherwise, you can choose not to normalize contact matrix by \"None\". "
	exit 1
    else
	echo "No normalization method is implemented."
    fi
fi

## *************************************************
## 4.2 Multi-reads can be reduced into uni-binPair
## if all the candiates alignment positions fall
## within on bin.
## *************************************************

if [ ! -d $dir/sorttmp ]; then
    mkdir -p $dir/sorttmp
fi

## take only the read ID and bin pairs information, then remove bin pair duplicates of one read pair
awk -v OFS="\t" 'BEGIN{q=0;c1=0;c2=0;b1=0;b2=0}(q!=$1 || c1!=$2 || c2!=$7 || b1!=$6 || b2!=$11){print $1, $2, $6, $7, $11;q=$1;c1=$2;c2=$7;b1=$6;b2=$11}' $validI.MULTI >$validI.MULTI.binPair


## Split multi-reads valid interaction bin pair files into unique bin pair and multiple bin pair file
awk '{print $1}' $validI.MULTI.binPair | sort -T $dir/sorttmp  | uniq -c | awk '$1 >1 {print $2}' >$validI.MULTI.binPair.multiList

rm -rf $dir/sorttmp
mkdir $dir/sorttmp


awk 'FNR==NR{a[$1];next};!($1 in a)' $validI.MULTI.binPair.multiList $validI.MULTI.binPair >$validI.MULTI.binPair.uni ##multi-reads reduced to uni-binPairs

awk 'FNR==NR{a[$1];next};($1 in a)' $validI.MULTI.binPair.multiList $validI.MULTI.binPair | sort -T $dir/sorttmp -k1,1V >$validI.MULTI.binPair.multi ## multi-reads fall into multi-binPairs

## uni bin pair count: Uni-binPair = Uni-reads binPair + Multi-reads reduced into uni-binPair
cat $validI.UNI.binPair $validI.MULTI.binPair.uni >$validI.binPair.uni
awk '{print $2, $3, $4, $5}' $validI.binPair.uni | sort -T $dir/sorttmp | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' | sort -k1,1V -k2,2n -k3,3V -k4,4n -T $dir/sorttmp >$validI.binPairCount.uni

rm -rf $dir/sorttmp
mkdir $dir/sorttmp


# Summary of removing duplicates
if [ "${#summaryFile}" -eq "0" ]; then
    summaryFile=$dir/rmDuplicates.summary
fi
echo -e "\n" >>$summaryFile
echo -e "\n" >>$summaryFile
echo -e "Step: Remove duplicates alignment pairs and binning" >>$summaryFile

validI_nodup_uni=$(cat $validI.UNI.binPair | wc -l)
validI_nodup_multi=$(cat $validI.MULTI | wc -l)
readCount_uniBin=$(cat $validI.binPair.uni | wc -l)
readCount_multiBin=$(cat $validI.MULTI.binPair.multiList | wc -l)


echo -e "Valid interaction pairs including uni-read pairs and all kinds of multi-mapping possible pairs:" >>$summaryFile
echo -e "  Valid interaction uni-mapping pairs noduplicates count:\t"$validI_nodup_uni >>$summaryFile
echo -e "  Valid interaction multi-mapping all possible pairs noduplicates count:\t"$validI_nodup_multi >>$summaryFile
echo -e "After binning, bin pair by read ID count:" >>$summaryFile
echo -e "  Valid read pair count with unique bin pair including uniquely determined multi-mapping read pairs:\t"$readCount_uniBin >>$summaryFile
echo -e "  Valid read pair count with multiple possible bin pairs:\t"$readCount_multiBin >>$summaryFile

## **********************************************
## 4.3 Uni-binPair contact matrix normalization.
## In preparation for step 5 prior building.
## **********************************************

# Exclude low count contacts
if [ "$minCount" -gt "1" ];then

    awk -v minC=$minCount -v OFS="\t" '$5>=minC {print $0}' $validI.binPairCount.uni | sort -k1,1V -k2,2n -k3,3V -k4,4n -T $dir/sorttmp  >$validI.binPairCount.uni.minCount$minCount

    if [ "$normMethod" != "None" ] &&  [ "$normMethod" != "NONE" ] &&  [ "$normMethod" != "none" ] &&  [ "$normMethod" != "" ]; then
	if [ "$normMethod" == "ICE" ] || [ "$normMethod" == "ice" ]; then

	    # ICE normalization with filtering low mappability regions
	    python3 $bin/ICE-with-sparseMatrix.py $validI.binPairCount.uni.minCount$minCount $mappFile l1 $validI.binPairCount.uni.ICEnorm $minMap $maxIter
	else
	    python3 $bin/KR_norm_mHiC.py -r $resolution -l $chromSizeFile -c $normChrom -tr $sparsePerc -f $validI.binPairCount.uni.minCount$minCount -o $dir
	fi
	## Uni bin marginal pair
	awk '{print $1, $2}' $validI.binPairCount.uni.${normMethod}norm >$validI.marginal1
	awk '{print $3, $4}' $validI.binPairCount.uni.${normMethod}norm >$validI.marginal2
    else
	## Uni bin marginal pair
	awk '{print $1, $2}' $validI.binPairCount.uni.minCount$minCount >$validI.marginal1
	awk '{print $3, $4}' $validI.binPairCount.uni.minCount$minCount >$validI.marginal2
    fi

else
    if [ "$normMethod" != "None" ] &&  [ "$normMethod" != "NONE" ] &&  [ "$normMethod" != "none" ] &&  [ "$normMethod" != "" ]; then

		if [ "$normMethod" == "ICE" ] || [ "$normMethod" == "ice" ]; then
		    # ICE normalization without filtering low mappability regions
		    # ICE-with-sparseMatrix.py just creates biases
		    echo "normalizing with ICE"
		    python3 $bin/ICE-with-sparseMatrix.py $validI.binPairCount.uni $mappFile l1 $validI.binPairCount.uni.ICEnorm $minMap $maxIter
		else
			echo "normalizing with KR"
		   python3 $bin/KR_norm_mHiC.py -r $resolution -l $chromSizeFile -c $normChrom -tr $sparsePerc -f $validI.binPairCount.uni -o $dir
		fi
		## Uni bin marginal pair
		awk '{print $1, $2}' $validI.binPairCount.uni.${normMethod}norm >$validI.marginal1
		awk '{print $3, $4}' $validI.binPairCount.uni.${normMethod}norm >$validI.marginal2
    else
		## Uni bin marginal pair
		echo "not normalizing LOL"
		awk '{print $1, $2}' $validI.binPairCount.uni >$validI.marginal1
		awk '{print $3, $4}' $validI.binPairCount.uni >$validI.marginal2
    fi
fi

cat $validI.marginal1 $validI.marginal2 | sort -k1,1V -k2,2n | awk '!a[$0]++' >$validI.binPair.Marginal
rm -rf $validI.marginal*
rm -rf $dir/sorttmp
