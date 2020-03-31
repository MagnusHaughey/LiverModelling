

import numpy as np 
import sys
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg




# First, obtain total number of cells, maxsize, from filename 
inFile = sys.argv[1]

maxSize = [string for string in inFile.split("/") if ("maxSize" in string)]
maxSize = [string for fragment in maxSize for string in fragment.split("_") if ("maxSize" in string)]
maxSize = int(maxSize[0].split("=")[1])


# Compute radius using maxSize as is done in the simulation script
radius_double = np.sqrt(maxSize/np.pi)
radius = int(2.0 * radius_double)



# Load data and re-scale x & y coordinates to radial distance 
x, y, stemCell, numMutations = np.loadtxt(inFile , unpack=True , delimiter=",")

r = [np.sqrt(((x-radius)**2) + ((y-radius)**2)) for x, y in zip(x, y)]



# Bin the data and compute mean values
xMin = 0
xMax = np.max(r)

yMin = 0
yMax = np.max(numMutations)

xBins = np.linspace(xMin, xMax, num=20)
yBins = np.linspace(yMin , yMax , num=500000)


hist, xedges, yedges = np.histogram2d(r , numMutations , (xBins, yBins))

hist = hist.T




f = open(sys.argv[1] + ".mutationsVSradialDistance.dat" , 'w')

toPlot_x , toPlot_weightedMean , toPlot_standardError = [] , [] , []


for LHboundary in range(len(xBins)-1):
	print("Processing horizontal bin {} of {}...".format( LHboundary , len(xBins)-1))
	binCentre_x = 0.5*(xBins[LHboundary] + xBins[LHboundary+1])

	# Compute weighted mean in y direction
	weightedMean = 0.0
	for bottomboundary in range(len(yBins)-1):
		binCentre_y = 0.5*(yBins[bottomboundary] + yBins[bottomboundary+1])
		weightedMean += binCentre_y * hist[bottomboundary][LHboundary]
		#print("{} {}".format( binCentre_y , hist[bottomboundary][LHboundary]) )

	if (np.sum(hist[:,LHboundary]) == 0):
		weightedMean = 0.0
	else:
		weightedMean /= np.sum(hist[:,LHboundary])


	# Compute standard error on mean 
	varSum = 0.0
	for bottomboundary in range(len(yBins)-1):
		binCentre_y = 0.5*(yBins[bottomboundary] + yBins[bottomboundary+1])
		varSum += hist[bottomboundary][LHboundary] * ((binCentre_y - weightedMean)**2)

	if (np.sum(hist[:,LHboundary]) - 1 == 0):
		varSum = 0.0
	else:
		varSum /= np.sum(hist[:,LHboundary]) - 1
		standardDeviation = np.sqrt(varSum)
		standardError = standardDeviation / np.sqrt(maxSize)


	f.write("{} {} {}\n".format(binCentre_x , weightedMean , standardError ))

	toPlot_x.append(binCentre_x)
	toPlot_weightedMean.append(weightedMean)
	toPlot_standardError.append(standardError)



f.close()




'''
toPlot_x , toPlot_weightedMean , toPlot_standardError = np.loadtxt(sys.argv[1] , unpack=True)



#plt.plot( toPlot_x , toPlot_weightedMean , c = 'Navy' , alpha = 0.2 , zorder = -10)
plt.errorbar( toPlot_x , toPlot_weightedMean , yerr=toPlot_standardError , fmt='' , ls = "None" , zorder = 0, capsize=4)
plt.scatter( toPlot_x , toPlot_weightedMean , c = 'Navy' , s = 12 , zorder = 10)
#plt.plot( toPlot_x , [val + err for val, err in zip(toPlot_weightedMean , toPlot_standardError)] , c = 'Navy' , alpha = 0.2 , zorder = -10)
#plt.plot( toPlot_x , [val - err for val, err in zip(toPlot_weightedMean , toPlot_standardError)] , c = 'Navy' , alpha = 0.2 , zorder = -10)


plt.tick_params(labelsize=16)
plt.xlabel(r"Radial distance, $r$" , fontsize=16)
plt.ylabel(r"Number of mutations" , fontsize=16)

plt.tight_layout()
plt.show()
'''


'''
# Image plotting

fig, (ax1 , ax2) = plt.subplots(1, 2 , figsize=(8,4) , gridspec_kw={'width_ratios': [2.5, 3]})


spatialDataImage = mpimg.imread(sys.argv[1] + '.png')
imgplot = ax1.imshow(spatialDataImage)
ax1.axis('off')

ax2.scatter( toPlot_x , toPlot_weightedMean , c = 'Navy' , s = 8)
ax2.plot( toPlot_x , toPlot_weightedMean , c = 'Navy' , alpha = 0.2 , zorder = -10)

ax2.tick_params(labelsize=16)
ax2.set_xlabel(r"Radial distance, $r$" , fontsize=16)
ax2.set_ylabel(r"Number of mutations" , fontsize=16)

# Manually set axis ranges (for continuity between different frames)
#ax2.set_xlim([0,165])
#ax2.set_ylim([0,55])

plt.tight_layout()

plt.subplots_adjust(wspace = 0.5)

#plt.show()
out_fig = sys.argv[1] + ".png.withScatterPlot.png"
plt.savefig(out_fig , dpi=500 , format='png')
'''




















