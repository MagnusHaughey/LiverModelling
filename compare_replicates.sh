

#/bin/bash 

for path in $1*A_summary.dat.corrected.dat;
do
	#echo $path

	file=$(echo $path | rev | cut -d "/" -f 1 | rev)
	#file=$(cut -d "/" -f 3 <<< $path)

	#echo $file

	IFS='_'
	read -ra fragments <<< "$file"

	IFS=' '
	primer=$(echo ${fragments[5]} | cut -c1-2)

	filepath=
	basic_filename=$1${fragments[0]}'_'${fragments[1]}'_'${fragments[2]}'_'${fragments[3]}'_'${fragments[4]}'_'$primer

	#echo $basic_filename
	

	#=========================	For somatic + germline data
	#replicateA=$basic_filename'A_summary.dat'
	#replicateB=$basic_filename'B_summary.dat'
	#outfile_name=$basic_filename'_shared_variants.dat'

	#replicateA=$basic_filename'A_summary.dat.corrected.dat'
	#replicateB=$basic_filename'B_summary.dat.corrected.dat'
	#outfile_name=$basic_filename'_shared_variants.dat.corrected.dat'



	#=========================	For germline removed data
	#replicateA=$basic_filename'A_summary.dat.noGermline.dat'
	#replicateB=$basic_filename'B_summary.dat.noGermline.dat'
	#outfile_name=$basic_filename'_shared_variants.dat.noGermline.dat'

	#replicateA=$basic_filename'A_summary.dat.noGermline.dat.corrected.dat'
	#replicateB=$basic_filename'B_summary.dat.noGermline.dat.corrected.dat'
	#outfile_name=$basic_filename'_shared_variants.dat.noGermline.dat.corrected.dat'


	if [[ $basic_filename != *"Stroma"* ]] && [[ ! -f $outfile_name ]] && [ -f "$replicateA" ] && [ -f "$replicateB" ]; then

		printf "Processing '$basic_filename'"

		python3 /Users/haughe01/Documents/LiverModelling/Year1/sequencing_analysis/scripts/SNV_checkReplicates.py \
				--replicateA $replicateA \
				--replicateB $replicateB \
				--output $basic_filename

	fi
	printf "\n"

done