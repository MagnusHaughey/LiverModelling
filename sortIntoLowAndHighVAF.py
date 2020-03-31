

import numpy as np 
import sys


#samplename = sys.argv[1].split("/")[-1].replace("GC_EC_8466_" , "").replace("_rep1_noGermline.dat" , "").replace("_rep1.dat" , "")
samplename = sys.argv[1].split("/")[-1].replace("GC_EC_8466_" , "").replace("_rep1_noGermline.dat.corrected.dat" , "").replace("_rep1.dat.corrected.dat" , "")

with open(sys.argv[1]) as f:
	ncols = len(f.readline().split('\t')) - 1


inData = np.genfromtxt(sys.argv[1] , dtype=float , skip_header=1 , delimiter="\t", usecols=range(1 , ncols))


inData = [['' if np.isnan(val) else val for val in row] for row in inData]


numClonal = 0
numLowVAFsubclonal = 0
numHighVAFsubclonal = 0

for row in inData:
	numVals = len([val for val in row if val != ''])
	if (numVals == (ncols - 1)):
		numClonal += 1
	else:
		numLowVAFsubclonal += len([val for val in row if val != '' and val < 0.1])
		numHighVAFsubclonal += len([val for val in row if val != '' and val >= 0.1])

print("{} {} {} {}".format(samplename , numClonal , numLowVAFsubclonal , numHighVAFsubclonal))

