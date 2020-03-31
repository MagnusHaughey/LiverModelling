

inFile=$1		# Data file. Should be file ending in ...bothPositionAndNucleotide.dat
outFile=$1.corrected.dat

refGenome="/Users/haughe01/Documents/LiverModelling/Year1/sequencing_analysis/scripts/NC_012920.1.fasta"


#========== For public & private mutations
#cat $inFile | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1 | sed 's/.\{2\}$//'

# Find sequence context
#for SNV in $(cat $inFile | cut -d" " -f1)
sed 1d $inFile | while read -r line
#while read -r line
do

	SNV=$(echo $line | cut -d" " -f1)
	
	pos=$(echo $SNV | sed 's/.\{2\}$//')
	base=${SNV: -2}

	printf "gi|251831106|ref|NC_012920.1|\t$((pos-1))\t$pos" > test.bed

	refBase=$(bedtools getfasta -fi $refGenome -bed test.bed | tail -1)
	#if [ $refBase == ${base:0:1} ]
	#then
	#	echo $line >> $outFile
	#fi

	if [ $refBase != ${base:0:1} ] && [ $refBase == ${base:1:1} ]
	then
		#echo $line
		#echo $pos
		echo $line | cut -d" " -f2-
		#while read -r freq; do
		#	echo $freq
		#done < <(echo $line | cut -d" " -f2-)
	fi

#done < $inFile
done



