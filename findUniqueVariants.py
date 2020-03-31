

import numpy as np 
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Read data
f = open(sys.argv[1] , 'r')

# Read lines of datafile
lines = f.readlines()

# Get rid of first line
lines = lines[1:]

# Strip each line of first and last entry
lines = [line.split('\t')[1:-1] for line in lines]

# Remove public variants
lines = [line for line in lines if not all([not element.isspace() for element in line])]

# Create binary matrix where 1 represents the presence of an SNV at particular point
binary_data = [['1' if not element.isspace() else '0' for element in line] for line in lines]


total_mutations = 0
for i in range(len(lines)):
	#print("{} --> {}".format(lines[i] , binary_data[i]))
	total_mutations += len([mut for mut in binary_data[i] if mut == '1'])


# Find all unique and non-unique variants for each cut
unique = np.zeros(len(lines[0]))
numPrivate = np.zeros(len(lines[0]))
for column in range(len(lines[0])):	# Loop over columns 
	for line in binary_data:	# Loop over rows

		variant_is_unique = True
		if (line[column] == '1'):

			for i in [-1,1]:	# If variant at line[column], then check neighbouring cuts

				#if (column - 1 >= 0) and (line[column - 1] == '0') and (column + 1 < len(lines[0])) and (line[column + 1] == '0')

				if (column + i < 0) or (column + i >= len(lines[0])): continue
				
				if (line[column + i] == '1'):
					variant_is_unique = False

			if variant_is_unique:
				unique[column] += 1

	# Express unique variants per cut as a fraction of total private variants
	numPrivate[column] = np.sum([int(line[column]) for line in binary_data])
	#print("C{} -> fraction private SNVs which are unique = {}".format( column+1 , unique[column]/numPrivate))





# Plot data
data_sample = sys.argv[1].split("_")[-(2+("noGermline" in sys.argv[1]))]
sample_names = ["{}C{}".format(data_sample , cut) for cut in range(1 , len(lines[0]) + 1)]

# For P4L3, remove the C2 bar (since only one amplicon was successful)
unique[1] = 0

plt.bar(np.linspace(1 , len(lines[0]) , len(lines[0])) , [unique/total*100 if total > 0 else 0 for unique , total in zip(unique , numPrivate)] , width=0.6 , facecolor=(1.0,0.656,0.388) , zorder = 10 , edgecolor='Black' , tick_label=sample_names)
#plt.bar(np.linspace(1 , len(lines[0]) , len(lines[0])) , [unique/total*100 if total > 0 else 0 for unique , total in zip(unique , numPrivate)] , width=0.6 , facecolor=(1.0,0.656,0.388) , zorder = 10 , edgecolor='Black' , tick_label=["Pa{}".format(sample[2:]) for sample in sample_names])

#plt.show()
plt.ylim([0,105])



# Turn on/off y label and ticks
#plt.ylabel('% unique private variants', fontsize=20)
ax = plt.gca()
plt.grid(axis='y' , zorder = -10 , alpha = 0.4)
ax.tick_params(axis='y', colors=(0,0,0,0))


#plt.grid(axis='y' , zorder = -10 , alpha = 0.4)

#plt.gca().tick_params(axis='both' , labelsize = 15)
plt.gca().tick_params(labelsize=20)
plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45 , ha="right")
plt.gca().set_yticklabels(['{:.0f}%'.format(y_tick) for y_tick in plt.gca().get_yticks()])

plt.tight_layout()

out_fig = sys.argv[1] + ".uniquePrivateVariants.png"
plt.savefig(out_fig , dpi=500 , format='png')








