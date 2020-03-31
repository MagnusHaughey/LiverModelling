

import numpy as np 
import sys
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr
from pylab import text


# Read data
clonal , lowVAFsubclonal , highVAFsubclonal = np.loadtxt(sys.argv[1] , unpack=True , delimiter = " " , usecols=(1,2,3))
clonal_noGermline , lowVAFsubclonal_noGermline , highVAFsubclonal_noGermline = np.loadtxt(sys.argv[2] , unpack=True , delimiter = " " , usecols=(1,2,3))

# Read in sample names
samplenames = np.genfromtxt(sys.argv[1] , unpack=True , usecols=(0) , dtype = str)
samplenames_noGermline = np.genfromtxt(sys.argv[2] , unpack=True , usecols=(0) , dtype = str)






#======================# Plot correlations

fig, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2, 3 , figsize=(18,9))
#fig.subplots_adjust(hspace=0.5, wspace=0.5)
ax1.set_xlim([-3,175])
ax2.set_xlim([-3,40])
ax3.set_xlim([-3,175])
ax4.set_xlim([-3,175])
ax5.set_xlim([-3,40])
ax6.set_xlim([-3,175])

ax1.set_ylim([-3,40])
ax2.set_ylim([-3,40])
ax3.set_ylim([-3,40])
ax4.set_ylim([-3,40])
ax5.set_ylim([-3,40])
ax6.set_ylim([-3,40])


#=== clonal vs low VAF subclonal

# Plot sample names by data points
for i in range(len(samplenames)):
	ax1.annotate(samplenames[i], (lowVAFsubclonal[i] + 1, clonal[i]) , fontsize=4)

# Plot correlations
ax1.scatter( lowVAFsubclonal , clonal , c = 'Black' , s = 4)
ax1.set_xlabel('# low VAF private SNVs', fontsize=16)
ax1.set_ylabel('# clonal SNVs', fontsize=16)

corr_1, _ = pearsonr(lowVAFsubclonal , clonal)


ax1.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_1 , _) , (105,33) , fontsize=13)

line1 = np.polyfit(lowVAFsubclonal , clonal, 1)



#print(line1[0])
#print(line1[1])
plotLine1 = [line1[0]*val + line1[1] for val in range(int(np.min(lowVAFsubclonal)) , int(np.max(lowVAFsubclonal)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax1.plot(range(int(np.min(lowVAFsubclonal)) , int(np.max(lowVAFsubclonal))) , plotLine1 , zorder = -10 , alpha = 0.5 , linestyle = 'dashdot')



#=== clonal vs high VAF subclonal

# Plot sample names by data points
for i in range(len(samplenames)):
	ax2.annotate(samplenames[i], (highVAFsubclonal[i] + 1, clonal[i]) , fontsize=4)

# Plot correlations
ax2.scatter( highVAFsubclonal , clonal , c = 'Black' , s = 4)
ax2.set_xlabel('# high VAF private SNVs', fontsize=16)
ax2.set_ylabel('# clonal SNVs', fontsize=16)

corr_2, _ = pearsonr(highVAFsubclonal , clonal)

ax2.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_2 , _) , (22,33) , fontsize=13)

line2 = np.polyfit(highVAFsubclonal , clonal, 1)

#print(line1[0])
#print(line1[1])
plotLine2 = [line2[0]*val + line2[1] for val in range(int(np.min(highVAFsubclonal)) , int(np.max(highVAFsubclonal)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax2.plot(range(int(np.min(highVAFsubclonal)) , int(np.max(highVAFsubclonal))) , plotLine2 , zorder = -10 , alpha = 0.5 , linestyle = 'dashdot')



#=== low VAF subclonal vs high VAF subclonal

# Plot sample names by data points
for i in range(len(samplenames)):
	ax3.annotate(samplenames[i], (lowVAFsubclonal[i] + 1, highVAFsubclonal[i]) , fontsize=4)

# Plot correlations
ax3.scatter( lowVAFsubclonal , highVAFsubclonal , c = 'Black' , s = 4)
ax3.set_xlabel('# low VAF private SNVs', fontsize=16)
ax3.set_ylabel('# high VAF private SNVs', fontsize=16)

corr_3, _ = pearsonr(lowVAFsubclonal , highVAFsubclonal)

ax3.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_3 , _) , (105,33) , fontsize=13)

line3 = np.polyfit(lowVAFsubclonal , highVAFsubclonal, 1)

#print(line1[0])
#print(line1[1])
plotLine3 = [line3[0]*val + line3[1] for val in range(int(np.min(lowVAFsubclonal)) , int(np.max(lowVAFsubclonal)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax3.plot(range(int(np.min(lowVAFsubclonal)) , int(np.max(lowVAFsubclonal))) , plotLine3 , zorder = -10 , alpha = 0.5 , linestyle = 'dashdot')







#=============================================================================








# Plot sample names by data points
for i in range(len(samplenames_noGermline)):
	ax4.annotate(samplenames_noGermline[i], (lowVAFsubclonal_noGermline[i] + 1, clonal_noGermline[i]) , fontsize=4)

# Plot correlations
ax4.scatter( lowVAFsubclonal_noGermline , clonal_noGermline , c = 'Black' , s = 4)
ax4.set_xlabel('# low VAF private SNVs', fontsize=16)
ax4.set_ylabel('# clonal SNVs', fontsize=16)

corr_1_noGermline, _ = pearsonr(lowVAFsubclonal_noGermline , clonal_noGermline)

ax4.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_1_noGermline , _) , (105,33) , fontsize=13)

line1_noGermline = np.polyfit(lowVAFsubclonal_noGermline , clonal_noGermline, 1)



#print(line1[0])
#print(line1[1])
plotLine1_noGermline = [line1_noGermline[0]*val + line1_noGermline[1] for val in range(int(np.min(lowVAFsubclonal_noGermline)) , int(np.max(lowVAFsubclonal_noGermline)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax4.plot(range(int(np.min(lowVAFsubclonal_noGermline)) , int(np.max(lowVAFsubclonal_noGermline))) , plotLine1_noGermline , zorder = -10 , alpha = 0.5 , color = "Red" , linestyle = 'dashdot')



#=== clonal vs high VAF subclonal

# Plot sample names by data points
for i in range(len(samplenames_noGermline)):
	ax5.annotate(samplenames_noGermline[i], (highVAFsubclonal_noGermline[i] + 1, clonal_noGermline[i]) , fontsize=4)

# Plot correlations
ax5.scatter( highVAFsubclonal_noGermline , clonal_noGermline , c = 'Black' , s = 4)
ax5.set_xlabel('# high VAF private SNVs', fontsize=16)
ax5.set_ylabel('# clonal SNVs', fontsize=16)

corr_2_noGermline, _ = pearsonr(highVAFsubclonal_noGermline, clonal_noGermline)

ax5.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_2_noGermline , _) , (22,33) , fontsize=13)

line2_noGermline = np.polyfit(highVAFsubclonal_noGermline , clonal_noGermline, 1)

#print(line1[0])
#print(line1[1])
plotLine2_noGermline = [line2_noGermline[0]*val + line2_noGermline[1] for val in range(int(np.min(highVAFsubclonal_noGermline)) , int(np.max(highVAFsubclonal_noGermline)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax5.plot(range(int(np.min(highVAFsubclonal_noGermline)) , int(np.max(highVAFsubclonal_noGermline))) , plotLine2_noGermline , zorder = -10 , alpha = 0.5 , color = "Red" , linestyle = 'dashdot')



#=== low VAF subclonal vs high VAF subclonal

# Plot sample names by data points
for i in range(len(samplenames_noGermline)):
	ax6.annotate(samplenames_noGermline[i], (lowVAFsubclonal_noGermline[i] + 1, highVAFsubclonal_noGermline[i]) , fontsize=4)

# Plot correlations
ax6.scatter( lowVAFsubclonal_noGermline , highVAFsubclonal_noGermline , c = 'Black' , s = 4)
ax6.set_xlabel('# low VAF private SNVs', fontsize=16)
ax6.set_ylabel('# high VAF private SNVs', fontsize=16)

corr_3_noGermline, _ = pearsonr(lowVAFsubclonal_noGermline , highVAFsubclonal_noGermline)

ax6.annotate("Pearson R = {0:.3f}\np-value = {1:.3f}".format(corr_3_noGermline , _) , (105,33) , fontsize=13)

line3_noGermline = np.polyfit(lowVAFsubclonal_noGermline , highVAFsubclonal_noGermline, 1)

#print(line1[0])
#print(line1[1])
plotLine3_noGermline = [line3_noGermline[0]*val + line3_noGermline[1] for val in range(int(np.min(lowVAFsubclonal_noGermline)) , int(np.max(lowVAFsubclonal_noGermline)))]

#plt.scatter(lowVAFsubclonal , plotLine1)
ax6.plot(range(int(np.min(lowVAFsubclonal_noGermline)) , int(np.max(lowVAFsubclonal_noGermline))) , plotLine3_noGermline , zorder = -10 , alpha = 0.5 , color = "Red" , linestyle = 'dashdot')


fig.tight_layout()
#plt.show()


out_fig = "correlations.png.corrected.png"
plt.savefig(out_fig , dpi=500 , format='png')



















