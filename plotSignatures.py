

import numpy as np
import sys
import matplotlib
matplotlib.use('ps')
from matplotlib import rc
import matplotlib.ticker as mtick
from matplotlib.ticker import FuncFormatter
import subprocess

rc('text',usetex=True)
rc('text.latex', preamble='\\usepackage{color} \\usepackage{xcolor} \
	\\definecolor{myOrange}{rgb}{1,0.5,0} \
	\\definecolor{myGreen}{rgb}{0.4864864864864865,0.8063063063063063,0.3738738738738739} \
	\\definecolor{myLightSeaGreen}{rgb}{0.08108108108108109,0.7747747747747747,0.6846846846846847} \
	\\definecolor{myNavy}{rgb}{0.11261261261261261,0.02252252252252252,0.7297297297297297} \
	\\definecolor{myPurple}{rgb}{0.7747747747747747,0.0945945945945946,0.9099099099099099}')

import matplotlib.pyplot as plt 

# Read datafile 
context , var = np.genfromtxt(sys.argv[1] , dtype = 'str' , unpack=True)

neighbouringBases = ['A','C','G','T']
mutations = ['CA' , 'CG' , 'CT' , 'TA' , 'TC' , 'TG']



# Sort mutations into bins 
counts = np.zeros((6,4,4))

for x in range(len(mutations)):
	for i in range(len(neighbouringBases)):
		for j in range(len(neighbouringBases)):

			for k in range(len(context)):
				if (context[k][1] == mutations[x][0]) and (var[k] == mutations[x][1]) and (context[k][0] == neighbouringBases[i]) and (context[k][-1] == neighbouringBases[j]):
					counts[x][i][j] += 1





labels = []
colours = []
allCounts = []
colourList = ['red' , (1,0.5,0) , (0.4864864864864865,0.8063063063063063,0.3738738738738739) , (0.08108108108108109,0.7747747747747747,0.6846846846846847) , (0.11261261261261261,0.02252252252252252,0.7297297297297297) , (0.7747747747747747,0.0945945945945946,0.9099099099099099)]
colourListLatex = ['red' , 'myOrange' , 'myGreen' , 'myLightSeaGreen' , 'blue' , 'myPurple']

for x in range(len(mutations)):

	for i in range(len(neighbouringBases)):
		for j in range(len(neighbouringBases)):

			labels.append("{}{}{}".format(neighbouringBases[i] , r'\textcolor{{{0}}}{{{1}}}'.format(colourListLatex[x] , mutations[x][0]) , neighbouringBases[j]))
			allCounts.append(counts[x][i][j])
			colours.append(colourList[x])


f = plt.figure(figsize=(8,2))
fig1 = f.add_subplot(111)
barlist = fig1.bar(np.linspace(1 , len(labels) , len(labels)) , [count/np.sum(allCounts)*100.0 for count in allCounts] , width=0.4 , color = colours , facecolor='Cyan' , edgecolor='none' , tick_label=labels)
for x in range(len(barlist)):
	barlist[x].set_color(colours[x])
plt.setp(fig1.xaxis.get_majorticklabels(), fontsize=5, rotation=90 )

plt.ylabel("Probability" , fontsize = 17 , weight='bold')
plt.yticks([5,10], fontsize=16)
plt.ylim([0,11])

f.set_tight_layout(True)

ax = plt.gca()
plt.gca().set_yticklabels(['{:.0f}\%'.format(x) for x in plt.gca().get_yticks()]) 
#ax.yaxis.set_major_formatter(mtick.PercentFormatter())


for i in range(len(mutations)):
	x = i*16 + 4
	y = 11.5
	plt.text(x, y, r"{} $\rightarrow$ {}".format(mutations[i][0] , mutations[i][1]), weight='extra bold', color=colourList[i], fontsize=12)

plt.text(90,9.5,r"$n = {}$".format(len(var)) , weight='extra bold', fontsize=12 )


out_fig = sys.argv[1] + ".ps"
plt.savefig(out_fig ,  dpi=500 )


# Convert ps to png using imageMagick
subprocess.call(["convert" , "-density" , "500" , "-geometry" , "100%", out_fig , sys.argv[1] + ".png"])
subprocess.call(["rm" , out_fig])









