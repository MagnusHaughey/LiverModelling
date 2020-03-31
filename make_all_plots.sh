

allFiles=$1'/*.csv'
allImages=$1'/%03d.csv.png.withScatterPlot.png'
outPath=$1'/out_withScatterPlot.mp4'

# Find maximum number of mutations in all cells in final timepoint
finalTimeData=$(ls -lrt $allFiles | tail -1 | rev | cut -d" " -f 1 | rev)
maxMutations=$(cat $finalTimeData | cut -d"," -f 4 | sort -nrk1,1 | head -1)


# Plot all spatial data
for file in $allFiles
do

	echo $file
	gnuplot -c ./gnuplot_spatialData.gp $file $maxMutations
	python3 plotNumberMutationsPerRadialDistance.py $file

done


# Combine all images into movie 
ffmpeg -framerate 10 -i $allImages -c:v libx264 -r 30 -pix_fmt yuv420p $outPath