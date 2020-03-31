
#/bin/bash 

for file in $1*rawDeepSNV.dat*
do
	trimmed_file=${file/_rawDeepSNV.dat/}

	rawData=$trimmed_file'_rawDeepSNV.dat'
	SNV_summary=$trimmed_file'_summary.dat'

	num_sig_variants_rep_A=$(grep "gi|251831106|ref|NC_012920.1|"  $rawData | wc -l)
	grep "gi|251831106|ref|NC_012920.1|"  $rawData > $(printf $rawData'.temp_one')
	grep "n.tst.fw" -A $num_sig_variants_rep_A $rawData > $(printf $rawData'.temp_two')
	grep "cov.ctrl.fw" -A $num_sig_variants_rep_A $rawData > $(printf $rawData'.temp_three')

	python3 /Users/haughe01/Documents/LiverModelling/Year1/sequencing_analysis/scripts/summariseSNVs_rCRS.py --input_one $(printf $rawData'.temp_one') --input_two $(printf $rawData'.temp_two') --input_three $(printf $rawData'.temp_three') --output $SNV_summary 

	rm $(printf $rawData'.temp_one') $(printf $rawData'.temp_two') $(printf $rawData'.temp_three')


done
