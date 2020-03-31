
import sys

# Ensure python can see numpy etc. packages
sys.path.append("/data/BCI-EvoCa2/magnus/liver_mtDNA/mtDNApipeline/.venv/lib/python3.6/site-packages/")


import numpy as np
import matplotlib.pyplot as plt
from pylab import text
import argparse
from operator import itemgetter




# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('-O', help='File path for outputted .png', type=str)
parser.add_argument('--data_origin', help='Specify raw coverage data obtained from samtools or picard', type=str)
parser.add_argument('--primer', help='Specify if primer M1 or M2, default is both M1+M2 combined', type=str, default='')
args = parser.parse_args()


if not ( (args.data_origin == 'SAMTOOLS') or (args.data_origin == 'PICARD') ): 
	raise Exception('--data_origin argument should be either "SAMTOOLS" or "PICARD".')

if not ( (args.primer == 'M1') or (args.primer == 'M2') or (args.primer == '') ): 
	raise Exception('--primer argument should be either "M1", "M2", or not given.')


# Obtain patient, sample and primer variables from output 
params = args.O.split("/")[-1]
params = params.strip("GC_EC_8466_")
params = params.strip("_cumulative_coverage.png")
patient, sample, primer = params.split("_")



# Import data from stdin
coverage_raw_data = sys.stdin.read()
coverage_raw_data = coverage_raw_data.splitlines()
#print(coverage_raw_data[0])



# Read depth at each position from third column data into ' depth[position , depth] '
position_depth = []
depth = []
bins = []


if (args.data_origin == 'SAMTOOLS'):
	for i in range(len(coverage_raw_data)):
			line_split = coverage_raw_data[i].split("\t")
			position_depth.append([int(line_split[1]),int(line_split[2])])
			depth.append(int(line_split[2]))

else:
	for i in range(1,len(coverage_raw_data)):
			line_split = coverage_raw_data[i].split("\t")
			bins.append(int(line_split[0]))
			#position_depth.append([int(line_split[0]),int(line_split[1])])
			depth.append(int(line_split[1]))




# Generate binned data 
#max_bin = max(position_depth, key=lambda x: x[0])[0]		# maximum bin, usually set to 20,000

#print(depth)


if (args.data_origin == 'SAMTOOLS'):
	bins = np.linspace(0 , 40000 , 400)
	
	# if only one (i.e. either M1 or M2) specified, then only use depth data for the appropriate region
	if (args.primer == 'M1'):
		depth = [pair[1] for pair in position_depth if ((pair[0] < 1465) or (pair[0] > 9397))]

	if (args.primer == 'M2'):
		depth = [pair[1] for pair in position_depth if ((pair[0] > 920) and (pair[0] < 9796))]


	depth, bins = np.histogram(depth , bins=bins)
	bins = bins[:-1]


else: 
	# Picard counts coverage for whole genome. So cut first ~10 bins to avoid all the zero counts and counts from mis-aligned reads
	#bins = np.linspace(0 , max_bin , max_bin+1)
	depth = depth[10::]
	bins = bins[10::]
	depth = np.asarray(depth)



# Remove trailing zeros from depth list, and corresponding bin value
while depth[-1] == 0:
    depth = depth[:-1]
    bins = bins[:-1]


#print(len(depth))
#print(len(bins))


cumulative_coverage = np.zeros(len(bins))
mean_depth = 0.0

#for i in range(1,len(depth)):
for i in range(len(bins)):
	cumulative_coverage[i] = np.sum(depth[i:])
	mean_depth += depth[i]*bins[i]

mean_depth /= np.sum(depth)


print(np.sum(depth))


# Normlise
cumulative_coverage /= (np.max(cumulative_coverage)/100.0)


#============= Plotting stuff

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
plt.xticks([0,10000,20000,30000,40000])

plt.plot(bins , cumulative_coverage, color='blue')
axes = plt.gca()
axes.set_xlim([0,40000])
axes.set_ylim([0,100])
plt.xlabel('Depth of coverage', fontsize=16)
plt.ylabel('% of target', fontsize=16)

text(0.70, 0.90,'Mean depth = {}'.format(int(mean_depth)), fontsize=16, horizontalalignment='center', verticalalignment='center', transform = axes.transAxes)

plt.fill_between(bins , cumulative_coverage, y2=0)
plt.tight_layout()

# Export plot
out_fig = args.O
plt.savefig(out_fig , dpi=300 , format='png')




