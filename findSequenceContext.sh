

inFile=$1		# Data file. Should be file ending in ...bothPositionAndNucleotide.dat

refGenome="/Users/haughe01/Documents/LiverModelling/Year1/sequencing_analysis/scripts/NC_012920.1.fasta"

equal=0
inequal=0
inequal_both_to_ref_and_var=0


#========== For public & private mutations
#cat $inFile | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1 | sed 's/.\{2\}$//'

# Find sequence context
#for SNV in $(cat $inFile | cut -f1)
for SNV in $(cat $inFile)
do

	
	pos=$(echo $SNV | sed 's/.\{2\}$//')
	base=${SNV: -2}

	printf "gi|251831106|ref|NC_012920.1|\t$((pos-2))\t$((pos+1))" > test.bed

	#if [ $refBase != ${base:0:1} ]
	#then
	#	#printf "$pos -> ${base:0:1} | $refBase | ${base:1:1}"
	#	#printf "\n"
	#
	#	if [ $refBase == ${base:1:1} ]
	#	then
	#		inequal=$((inequal + 1))
	#	else 
	#		inequal_both_to_ref_and_var=$(( inequal_both_to_ref_and_var + 1 ))
	#	fi
	#else
	#	equal=$((equal + 1))
	#fi


	if [ ${base:0:1} == 'C' ]
	then
		if [ ${base:1:1} == 'A' ] || [ ${base:1:1} == 'G' ] || [ ${base:1:1} == 'T' ]
		then
			printf "$(bedtools getfasta -fi $refGenome -bed test.bed | tail -1) ${base:1:1}\n"
		fi
	fi
	if [ ${base:0:1} == 'T' ]
	then
		if [ ${base:1:1} == 'A' ] || [ ${base:1:1} == 'G' ] || [ ${base:1:1} == 'C' ]
		then
			printf "$(bedtools getfasta -fi $refGenome -bed test.bed | tail -1) ${base:1:1}\n"
		fi
	fi

done


#total_bases=$(cat $inFile | wc -l)


#printf "$inFile\t$total_bases\t$equal\t$inequal\t$inequal_both_to_ref_and_var\n"



