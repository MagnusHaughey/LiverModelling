
import numpy as np 
import sys
import os
#os.system('ls -l')

# Read in reference sequence
raw_chrM = np.genfromtxt("/Users/haughe01/Documents/LiverModelling/Year1/sequencing_analysis/scripts/NC_012920.1.fasta", dtype='str', skip_header=1)
chrM = [base for line in raw_chrM for base in line]

# Read mutations data
#mutations = np.genfromtxt(sys.argv[1] , dtype='str')

#positions = [int(mutation[:-2]) for mutation in mutations]

positions = np.loadtxt("./test.dat.txt", usecols=(2))
ref = np.genfromtxt("./test.dat.txt" , dtype='str' , usecols=(3))
var = np.genfromtxt("./test.dat.txt" , dtype='str' , usecols=(4))

#print([int(pos) for pos in positions])
#print(ref)

for i in range(len(positions)):
	print("{}{}{}".format(int(positions[i]) , ref[i] , var[i]))

#print(numEqual)
#print(numNotEqual)

exit(0)


numEqual = 0
numNotEqual = 0
for i in range(len(positions)):
	if (mutations[i][-2] == chrM[positions[i] - 1]):
		numEqual += 1
		print("{} {} {}".format(positions[i] , mutations[i][-2] , chrM[positions[i] - 1]))
	else:
		numNotEqual += 1


print(numEqual)
print(numNotEqual)