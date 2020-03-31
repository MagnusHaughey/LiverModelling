
import numpy as np 
import sys 
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

inData = np.genfromtxt(sys.argv[1] , 'str')

allSNVtypes = sorted(set(inData))

# Add this line in so that both PT plots have the same number of columns
#allSNVtypes.append("AC")
#allSNVtypes = sorted(set(allSNVtypes))

SNVrawNumbers = []

for SNV in allSNVtypes:
	SNVrawNumbers.append(len([var for var in inData if (var == SNV)]))
	#print("{}: {}".format( SNV , len([var for var in inData if (var == SNV)])))

totalNumberSNVs = np.sum(SNVrawNumbers)

sample_names = allSNVtypes
sample_names = [r'{}$\rightarrow${}'.format(var[0] , var[1]) for var in allSNVtypes]

barColours = []
for var in sample_names:
	if var[0] == 'A':
		barColours.append((1, 0.1, 0.1, 0.5))
	elif var[0] == 'C':
		barColours.append((1, 0.7, 0.2, 0.5))
	elif var[0] == 'G':
		barColours.append((0.2, 0.7, 0.3, 0.5))
	elif var[0] == 'T':
		barColours.append((0.2, 0.4, 0.7, 0.5))

barChart = plt.bar( range(len(SNVrawNumbers)) , [freq*100 for freq in SNVrawNumbers]/totalNumberSNVs , width=0.6 , color=barColours , edgecolor='Black' , tick_label=sample_names)
plt.xticks(fontsize=18 , rotation=90 , weight='bold')
plt.yticks(fontsize=15)
plt.minorticks_off()

plt.ylim([0,45])

# Remove y tics if necessary
#plt.yticks([], [])


axes = plt.gca()
axes.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
plt.tight_layout()

#plt.show()

out_fig = sys.argv[1] + ".png"
plt.savefig(out_fig , dpi=500 , format='png')



















































