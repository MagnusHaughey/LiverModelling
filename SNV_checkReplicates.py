
import numpy as np 
import sys
import argparse
import matplotlib.pyplot as plt 


# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('--replicateA', help='', type=str)
parser.add_argument('--replicateB', help='', type=str)
parser.add_argument('--output', help='', type=str)
args = parser.parse_args()

positionsA = np.genfromtxt(args.replicateA , unpack=True , usecols=(0) , dtype=str)
positionsB = np.genfromtxt(args.replicateB , unpack=True , usecols=(0) , dtype=str)

#positionsA = np.array([positionsA])
#positionsB = np.array([positionsB])

frequenciesA = np.loadtxt(args.replicateA , unpack=True , usecols=(1))
frequenciesB = np.loadtxt(args.replicateB , unpack=True , usecols=(1))




if (positionsA.size == 1):
	positionsA = [str(positionsA)]
	frequenciesA = [frequenciesA]

if (positionsB.size == 1):
	positionsB = [str(positionsB)]
	frequenciesB = [frequenciesB]



#print(positionsA.shape[0])
#print(positionsB.shape[0])
#print(type(frequenciesB))

shared_variants , positions_written , label , rep1 , rep2 = [] , [] , [] , [] , []

f = open(args.output + "_separate_replicate_frequencies.dat.corrected.dat" , 'w')
#for i in range(len(positionsA)):

for i in range(len(positionsA)):
	
	shared = False
	#for j in range(len(positionsB)):
	for j in range(len(positionsB)):

		if (positionsB[j] == positionsA[i]):
			shared = True

			label.append("{} {}>{}".format(positionsA[i][:-2] , positionsA[i][-2] , positionsA[i][-1]))
			rep1.append(frequenciesA[i])
			rep2.append(frequenciesB[j])

			f.write("{}{}{}\t{}\t{}\n".format(positionsA[i][:-2] , positionsA[i][-2] , positionsA[i][-1] , frequenciesA[i] , frequenciesB[j]))
			positions_written.append('{}{}{}'.format(positionsA[i][:-2] , positionsA[i][-2] , positionsA[i][-1]))

for i in range(len(positionsA)):
	if ('{}{}{}'.format(positionsA[i][:-2] , positionsA[i][-2] , positionsA[i][-1]) not in positions_written):
		f.write("{}{}{}\t{}\t\n".format(positionsA[i][:-2] , positionsA[i][-2] , positionsA[i][-1] , frequenciesA[i]))

for i in range(len(positionsB)):
	if ('{}{}{}'.format(positionsB[i][:-2] , positionsB[i][-2] , positionsB[i][-1]) not in positions_written):
		f.write("{}{}{}\t\t{}\n".format(positionsB[i][:-2] , positionsB[i][-2] , positionsB[i][-1] , frequenciesB[i]))

f.close()


plt.plot([0, 1], [0, 1], color = 'red', linewidth = 1.5 , zorder=-1)
plt.scatter(rep1 , rep2 , s=2 , color='black' , zorder=1)
axes = plt.gca()
axes.set_xlim([0.0001,1])
axes.set_ylim([0.0001,1])

for i in range(len(label)):
    plt.annotate(label[i], (rep1[i]+np.log(1-rep1[i]*0.4)*0.4, rep2[i]+np.log(1-rep2[i]*0.4)*0.4) , fontsize=4)

plt.yscale("log")
plt.xscale("log")

plt.xlabel("Frequency in repeat 1")
plt.ylabel("Frequency in repeat 2")

outfig = args.output + "_replicate_frequencies.png.corrected.png"
plt.savefig(outfig , dpi=500 , format='png')


'''
f = open(args.output + "_shared_variants.dat" , 'w')
for i in range(len(rep1)):
	f.write("{} {} {}\n".format(label[i] , rep1[i] , rep2[i]))

f.close()
'''



