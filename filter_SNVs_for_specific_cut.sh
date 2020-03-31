

inFile=$1		# Data file. Should be final processed SNV data (i.e. file ending in ...rep1.dat or ...rep2.dat)
cut=$(( $2 + 1 ))	# Cut number for which you want to print mutation positions
#fasta=$3




#========== For public & private mutations
#cat $inFile | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1 | sed 's/.\{2\}$//'

# Find sequence context
#for SNV in $(cat $inFile | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1)
#do
#	
#	printf "$(echo $SNV | sed 's/.\{2\}$//')${SNV: -2}\n"
#done


# Print positions of mutations (public & private) -> nucleotides only 
#for SNV in $(cat $inFile | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1)
#do
#	#echo ${SNV: -2}
#	echo $SNV
#done




#========== For private mutations only
numColumns=$(awk '{print NF}' $inFile | sort -nu | tail -n 1)

while read line
do
	numMutationsInLine=$(echo $line | awk '{print NF}' | sort -nu | tail -n 1)

	if (( numMutationsInLine == numColumns ))
	then

		lastPublicMutation=$line
		#echo $line

	fi
	#printf "$line ||| $numMutationsInLine \n"

done < $inFile

startFrom=$(grep -Fn "$lastPublicMutation" $inFile | cut -d":" -f1)
totalNumLines=$(cat $inFile | wc -l)


# Print positions of mutations (private only) without nucleotides
#cat $inFile | tail -n $((totalNumLines - startFrom)) | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1 | sed 's/.\{2\}$//'

# Print positions of mutations (private only) -> nucleotides only 
for SNV in $(cat $inFile | tail -n $((totalNumLines - startFrom)) | cut -d$'\t' -f1,$cut | grep "\." | cut -d$'\t' -f1)
do
	#echo ${SNV: -2}
	echo $SNV
done
