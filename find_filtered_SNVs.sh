

#/bin/bash 

inFileMetrics=$1'.METRICS.dat'

for variant in $(cat $inFileMetrics | cut -d" " -f 3)
do

	otherReplicate={$1/M1A/M1B}

	otherReplicateMetrics=$otherReplicate'.METRICS.dat'

	grep $variant $otherReplicate $otherReplicateMetrics

done