


import numpy as np 
import matplotlib.pyplot as plt
import glob
import sys
import os




def average_double_entries(pairs):


	# Identify duplicates
	duplicates = []
	for position in [pos[0] for pos in pairs]:
		if ([pos[0] for pos in pairs].count(position) > 1) and (position not in duplicates):
			duplicates.append(position)

	# Remove duplicates
	means = []
	for dup in duplicates:
		freqs = []
		for j in range(len(pairs)):
			if (pairs[j][0] == dup):
				freqs.append(float(pairs[j][1]))

		means.append(np.mean(freqs))


	# Remove all instances of duplicate entry, and then replace with single, averaged entry
	for j in range(len(duplicates)):
		pairs = [pair for pair in pairs if (pair[0] != duplicates[j])]
		pairs.append((duplicates[j] , str(means[j])))


	return pairs




parent_dir = sys.argv[1]

all_files = []
#for file in glob.glob(parent_dir + "*summary.dat"):		# should be specified up to the patient name, e.g. (/path/GC_EC_8466_193B)
#	all_files.append(file[:-16])


all_files = []
for file in glob.glob(parent_dir + "*summary.dat.noGermline.dat.corrected.dat"):		# should be specified up to the patient name, e.g. (/path/GC_EC_8466_193B)
	all_files.append(file.replace("_summary.dat.noGermline.dat.corrected.dat" , "")[:-4])


all_files = sorted(set(all_files))


split_filenames = [file.split("_") for file in all_files if ("Bulk" not in file)]
#split_filenames = [file.split("_") for file in all_files]
sample_names = [frag[-1] for frag in split_filenames]




all_repA , all_repB = [] , []
for i in range(len(all_files)):

	if ("Bulk" in all_files[i]): 
		continue

	missing_primer = ''
	counterpart = ''

	f = open("errors.dat" , 'a')
	# First, check for files
	for primer in ['M1A' , 'M1B' , 'M2A' , 'M2B']:
		#fileInfo = os.stat(all_files[i] + "_" + primer + "_summary.dat")
		fileInfo = os.stat(all_files[i] + "_" + primer + "_summary.dat.noGermline.dat.corrected.dat")
		if (fileInfo.st_size == 0):
			f.write("*Warning, file '" + all_files[i] + "_" + primer + "_summary.dat' is empty.\n")
			f.flush()

			missing_primer = primer
			
			# Determine the counterpart to the missing primer i.e. M1A -> M1B, M2A -> M2B and vice versa
			if primer in ['M1A' , 'M1B']:
				primer_pair = ['M1A' , 'M1B']
				primer_pair.remove(primer)
				counterpart = primer_pair[0]
			else:
				primer_pair = ['M2A' , 'M2B']
				primer_pair.remove(primer)
				counterpart = primer_pair[0]


			# If corresponding second replicate is also missing, then abort
			#fileInfo = os.stat(all_files[i] + "_" + counterpart + "_summary.dat")
			fileInfo = os.stat(all_files[i] + "_" + counterpart + "_summary.dat.noGermline.dat.corrected.dat")
			if (fileInfo.st_size == 0):
				print("Both `" + all_files[i] + "_" + primer + "_summary.dat` and `" + all_files[i] + "_" + counterpart + "_summary.dat` are missing. Exiting.")
				#exit(0)
			#print("Careful when processing primer " + counterpart)

	f.close()


	if (missing_primer != 'M1A'):
		#M1A_pos , M1A_freq = np.genfromtxt(all_files[i] + "_M1A_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M1A_pos , M1A_freq = np.genfromtxt(all_files[i] + "_M1A_summary.dat.noGermline.dat.corrected.dat" , unpack=True , usecols=(0,1) , dtype=str)
		if (isinstance(M1A_pos, str)): # If data file only contains one variant (single line) then results of np.genfromtxt will be a string, instead of a list of strings. Avoid errors by casting single value into a list 
			M1A_pos = [M1A_pos]
			M1A_freq = [M1A_freq]
		M1A_position_frequency = [ (pos , freq) for pos, freq in zip(M1A_pos , M1A_freq)]
	else:
		# if primer is missing, then fill its arrays with data for *other* replicate. This is to ensure that the program runs correctly
		M1A_position_frequency = []
		#M1A_pos , M1A_freq = np.genfromtxt(all_files[i] + "_M1B_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		#M1A_position_frequency = [ (pos , freq) for pos, freq in zip(M1B_pos , M1B_freq)]

	if (missing_primer != 'M1B'):
		#M1B_pos , M1B_freq = np.genfromtxt(all_files[i] + "_M1B_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M1B_pos , M1B_freq = np.genfromtxt(all_files[i] + "_M1B_summary.dat.noGermline.dat.corrected.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M1B_position_frequency = [ (pos , freq) for pos, freq in zip(M1B_pos , M1B_freq)]
	else:
		M1B_position_frequency = []
		#M1B_pos , M1B_freq = np.genfromtxt(all_files[i] + "_M1A_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		#M1B_position_frequency = [ (pos , freq) for pos, freq in zip(M1B_pos , M1B_freq)]

	if (missing_primer != 'M2A'):
		#M2A_pos , M2A_freq = np.genfromtxt(all_files[i] + "_M2A_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M2A_pos , M2A_freq = np.genfromtxt(all_files[i] + "_M2A_summary.dat.noGermline.dat.corrected.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M2A_position_frequency = [ (pos , freq) for pos, freq in zip(M2A_pos , M2A_freq)]
	else:
		M2A_position_frequency = []
		#M2A_pos , M2A_freq = np.genfromtxt(all_files[i] + "_M2B_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		#M2A_position_frequency = [ (pos , freq) for pos, freq in zip(M2A_pos , M2A_freq)]

	if (missing_primer != 'M2B'):
		#M2B_pos , M2B_freq = np.genfromtxt(all_files[i] + "_M2B_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M2B_pos , M2B_freq = np.genfromtxt(all_files[i] + "_M2B_summary.dat.noGermline.dat.corrected.dat" , unpack=True , usecols=(0,1) , dtype=str)
		M2B_position_frequency = [ (pos , freq) for pos, freq in zip(M2B_pos , M2B_freq)]
	else:
		M2B_position_frequency = []
		#M2B_pos , M2B_freq = np.genfromtxt(all_files[i] + "_M2A_summary.dat" , unpack=True , usecols=(0,1) , dtype=str)
		#M2B_position_frequency = [ (pos , freq) for pos, freq in zip(M2B_pos , M2B_freq)]


	# Join data from two primers for each replicate
	repA_position_frequency = M1A_position_frequency + M2A_position_frequency
	repB_position_frequency = M1B_position_frequency + M2B_position_frequency


	# Exclude any variants that are not called in both A and B repeats
	for position in [pos[0] for pos in repA_position_frequency]:

		if (position not in [pos[0] for pos in repB_position_frequency]):
			index = [pos[0] for pos in repA_position_frequency].index(position)
			del repA_position_frequency[index]

	for position in [pos[0] for pos in repB_position_frequency]:

		if (position not in [pos[0] for pos in repA_position_frequency]):
			index = [pos[0] for pos in repB_position_frequency].index(position)
			del repB_position_frequency[index]


	# Occasionally where the primers overlap both M1 and M2 call the same variant. If so, average their frequencies 
	if (len(set([pos[0] for pos in repA_position_frequency])) != len([pos[0] for pos in repA_position_frequency])):
		repA_position_frequency = average_double_entries(repA_position_frequency)
		

	if (len(set([pos[0] for pos in repB_position_frequency])) != len([pos[0] for pos in repB_position_frequency])):
		repB_position_frequency = average_double_entries(repB_position_frequency)
 



	all_repA.append(repA_position_frequency)
	all_repB.append(repB_position_frequency)


#print(sample_names)




#==================== Find shared (clonal) mutations
outfile_name = all_files[0][:all_files[0].rindex('_')]
#f = open(parent_dir + "_rep1.dat" , 'w')
f = open(parent_dir + "_rep1_noGermline.dat.corrected.dat" , 'w')

f.write("\t")
for i in range(len(sample_names)):
	f.write("{:12.10s}\t".format(sample_names[i]))
f.write("\n")



#print("\n*** Shared mutations\n\n")
summed_frequencies = np.zeros(len(all_repA))
number_of_variants = np.zeros(len(all_repA))
shared_mutations = []
clonal_frequencies = []
for position in [pos[0] for pos in all_repA[0]]:

	shared = True
	for i in range(1,len(all_repA)):
		if (position not in [pos[0] for pos in all_repA[i]]):
			shared = False

	if shared:
		shared_mutations.append(position)
		freqs = []
		for i in range(len(all_repA)):
			freqs.append([pair[1] for pair in all_repA[i] if (pair[0] == position)])

		clonal_frequencies.append([position , freqs])

		#print("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
		#f.write("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}\n".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
		f.write("{}\t".format(position))
		for i in range(len(all_repA)):
			f.write("{:12.10s}\t".format(freqs[i][0]))
			summed_frequencies[i] += float(freqs[i][0])
			number_of_variants[i] += 1
		f.write("\n")


#==================== Private mutations

#print("\n*** Private mutations\n\n")

already_printed = []
for file in range(len(all_repA)):		# Loop over all files and print private mutations

	for position in [pos[0] for pos in all_repA[file]]:		

		if (position not in shared_mutations) and (position not in already_printed):
			already_printed.append(position)

			# Loop over all files and find frequency of this private mutation
			freqs = []
			for i in range(len(all_repA)):
				val = [pair[1] for pair in all_repA[i] if (pair[0] == position)]
				if (len(val) > 0):
					freqs.append([pair[1] for pair in all_repA[i] if (pair[0] == position)])
				else:
					freqs.append(" ")

			#print("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
			#f.write("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}\n".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
			f.write("{}\t".format(position))
			for i in range(len(all_repA)):
				f.write("{:12.10s}\t".format(freqs[i][0]))
				if (freqs[i][0] == ' '):
					summed_frequencies[i] += 0.0
				else:
					summed_frequencies[i] += float(freqs[i][0])
					number_of_variants[i] += 1
			f.write("\n")


f.close()


#f = plt.figure(figsize=(10,4))
#fig1 = f.add_subplot(121)
#fig2 = f.add_subplot(122)

f = plt.figure(figsize=(6,5))
fig1 = f.add_subplot(111)

#print(number_of_variants)
#print(summed_frequencies)

# For P4L3, remove the C2 bar (since only one amplicon was successful)
#number_of_variants[1] = 0

#fig1.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , number_of_variants , width=0.6 , facecolor='Cyan' , edgecolor='Black' , tick_label=sample_names)
fig1.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , number_of_variants , width=0.6 , facecolor='Red' , alpha=0.5 , edgecolor='Black' , tick_label=["PbC{}".format(sample) for sample in sample_names])
fig1.set_ylabel('No. variants', fontsize=25)

#fig2.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , summed_frequencies , width=0.6 , facecolor='Cyan' , edgecolor='Black' , tick_label=sample_names)
#fig2.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , summed_frequencies , width=0.6 , facecolor='Red' , alpha=0.5 , edgecolor='Black' , tick_label=sample_names)
#fig2.set_ylabel('Summed frequencies', fontsize=16)

f.gca().tick_params(labelsize=20)
#f.gca().set_yticks([0,5,10]) 
f.subplots_adjust(wspace=0.4, bottom=0.2)

plt.setp(fig1.xaxis.get_majorticklabels(), rotation=45 , ha="right")
#plt.setp(fig2.xaxis.get_majorticklabels(), rotation=45 , ha="right")

plt.tight_layout()

#out_fig = parent_dir + "_rep1_mutational_burden.png"
out_fig = parent_dir + "_rep1_mutational_burden_noGermline.png.corrected.png"
f.savefig(out_fig , dpi=500 , format='png')

plt.clf()

if (len(clonal_frequencies) > 0):

	# Plot all clonal variants
	for clonal in clonal_frequencies:

		line, = plt.plot(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , alpha = 0.5 , zorder = -10 , label=clonal[0])
		plt.scatter(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , s =10 , zorder = 0)

		#print([float(item) for sublist in clonal[1] for item in sublist])


	# If there are both high and low VAF clonal variants, control vertical range if needs be. Otherwise have the y-range set to [0:1]
	axes = plt.gca()
	axes.set_ylim([0,1])

	plt.xticks(np.arange(1 , len(clonal_frequencies[0][1]) + 1) , ["C{}".format(sample) for sample in sample_names] , fontsize=10 , rotation = 45 , ha="right")
	#plt.xlabel('Zone', fontsize=20)

	plt.yticks(fontsize=10)
	plt.ylabel('VAF', fontsize=20)

	plt.legend()

	out_fig = parent_dir + "_rep1_clonal_freqiencies_noGermline.png.corrected.png"
	plt.savefig(out_fig , dpi=500 , format='png' , bbox_inches='tight')
	plt.clf()

	#=================================================================================

	# Plot high VAF clonal variants 
	nonEmptyPlot = False
	for clonal in clonal_frequencies:

		if np.min([float(item) for sublist in clonal[1] for item in sublist]) > 0.2:

			line, = plt.plot(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , alpha = 0.5 , zorder = -10 , label=clonal[0])
			plt.scatter(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , s =10 , zorder = 0)
			nonEmptyPlot = True

		#print([float(item) for sublist in clonal[1] for item in sublist])

	if (nonEmptyPlot == True):

		# If there are both high and low VAF clonal variants, control vertical range if needs be. Otherwise have the y-range set to [0:1]
		axes = plt.gca()
		axes.set_ylim([0.5,1])

		#plt.xticks(np.arange(1 , len(clonal_frequencies[0][1]) + 1) , sample_names , fontsize=10 , rotation = 45 , ha="right")
		plt.xticks(np.arange(1 , len(clonal_frequencies[0][1]) + 1) , ["C{}".format(sample) for sample in sample_names] , fontsize=10 , rotation = 0)
		#plt.xlabel('Zone', fontsize=20)

		plt.yticks(fontsize=10)
		plt.ylabel('VAF', fontsize=20)

		plt.legend()

		out_fig = parent_dir + "_rep1_clonal_freqiencies_highVAF_noGermline.png.corrected.png"
		plt.savefig(out_fig , dpi=500 , format='png' , bbox_inches='tight')
		plt.clf()

	#=================================================================================

	# Plot low VAF clonal variants 
	nonEmptyPlot = False
	for clonal in clonal_frequencies:


		if np.max([float(item) for sublist in clonal[1] for item in sublist]) < 0.2:

			#print([float(item) for sublist in clonal[1] for item in sublist])

			line, = plt.plot(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , alpha = 0.5 , zorder = -10 , label=clonal[0])
			plt.scatter(np.arange(1 , len(clonal[1]) + 1), [float(item) for sublist in clonal[1] for item in sublist] , s =10 , zorder = 0)
			nonEmptyPlot = True

		#print([float(item) for sublist in clonal[1] for item in sublist])

	if (nonEmptyPlot == True):



		# If there are both high and low VAF clonal variants, control vertical range if needs be. Otherwise have the y-range set to [0:1]
		axes = plt.gca()
		axes.set_ylim([0,0.2])

		#plt.xticks(np.arange(1 , len(clonal_frequencies[0][1]) + 1) , sample_names , fontsize=10 , rotation=45 , ha="right")
		plt.xticks(np.arange(1 , len(clonal_frequencies[0][1]) + 1) , ["C{}".format(sample) for sample in sample_names] , fontsize=10 , rotation=0)
		#plt.xlabel('Zone', fontsize=20)

		plt.yticks(fontsize=10)
		plt.ylabel('VAF', fontsize=20)

		plt.legend()

		out_fig = parent_dir + "_rep1_clonal_freqiencies_lowVAF_noGermline.png.corrected.png"
		plt.savefig(out_fig , dpi=500 , format='png' , bbox_inches='tight')
		plt.clf()

#=====================================================================================================

#==================== Find shared (clonal) mutations


summed_frequencies = np.zeros(len(all_repB))
number_of_variants = np.zeros(len(all_repB))

#g = open(parent_dir + "_rep2.dat" , 'w')
g = open(parent_dir + "_rep2_noGermline.dat.corrected.dat" , 'w')


g.write("\t")
for i in range(len(sample_names)):
	g.write("{:12.10s}\t".format(sample_names[i]))
g.write("\n")

#print("\n\nRep 2\n")
for position in [pos[0] for pos in all_repB[0]]:

	shared = True
	for i in range(1,len(all_repB)):
		if (position not in [pos[0] for pos in all_repB[i]]):
			shared = False

	if shared:
		freqs = []
		for i in range(len(all_repB)):
			freqs.append([pair[1] for pair in all_repB[i] if (pair[0] == position)])

		#print("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
		#g.write("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}\n".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
		g.write("{}\t".format(position))
		for i in range(len(all_repB)):
			g.write("{:12.10s}\t".format(freqs[i][0]))
			summed_frequencies[i] += float(freqs[i][0])
			number_of_variants[i] += 1
		g.write("\n")



#==================== Private mutations

#print("\n*** Private mutations\n\n")

#already_printed = []
for position in already_printed:

	freqs = []
	for i in range(len(all_repB)):
		val = [pair[1] for pair in all_repA[i] if (pair[0] == position)]
		if (len(val) > 0):
			freqs.append([pair[1] for pair in all_repB[i] if (pair[0] == position)])
		else:
			freqs.append(" ")

	#print("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))
	g.write("{}\t".format(position))
	for i in range(len(all_repB)):
		g.write("{:12.10s}\t".format(freqs[i][0]))
		if (freqs[i][0] == ' '):
			summed_frequencies[i] += 0.0
		else:
			summed_frequencies[i] += float(freqs[i][0])
			number_of_variants[i] += 1
	g.write("\n")
	#g.write("{}\t {:12.10s}\t {:12.10s}\t {:12.10s}\t {:12.10s}\n".format(position , freqs[0][0] , freqs[1][0] , freqs[2][0] , freqs[3][0]))



g.close()


# Plot summed mutation frequency data 

#f = plt.figure(figsize=(10,4))
#fig1 = f.add_subplot(121)
#fig2 = f.add_subplot(122)

f = plt.figure(figsize=(5,4))
fig1 = f.add_subplot(111)

#print(number_of_variants)
#print(summed_frequencies)

#fig1.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , number_of_variants , width=0.6 , facecolor='Cyan' , edgecolor='Black' , tick_label=sample_names)
fig1.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , number_of_variants , width=0.6 , facecolor='Red' , alpha=0.5 , edgecolor='Black' , tick_label=["C{}".format(sample) for sample in sample_names])
fig1.set_ylabel('No. variants', fontsize=16)

#fig2.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , summed_frequencies , width=0.6 , facecolor='Cyan' , edgecolor='Black' , tick_label=sample_names)
#fig2.bar(np.linspace(1 , len(all_repA) , len(all_repA)) , summed_frequencies , width=0.6 , facecolor='Red' , alpha=0.5 , edgecolor='Black' , tick_label=sample_names)
#fig2.set_ylabel('Summed frequencies', fontsize=16)

f.subplots_adjust(wspace=0.4, bottom=0.2)

#plt.setp(fig1.xaxis.get_majorticklabels(), rotation=45 , ha="right")
#plt.setp(fig2.xaxis.get_majorticklabels(), rotation=45 , ha="right")

plt.tight_layout()

#out_fig = parent_dir + "_rep2_mutational_burden.png"
out_fig = parent_dir + "_rep2_mutational_burden_noGermline.png.corrected.png"
f.savefig(out_fig , dpi=500 , format='png')






