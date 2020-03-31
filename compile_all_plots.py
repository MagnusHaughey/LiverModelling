
import numpy as np 
import glob
import subprocess
import sys

parent_dir = sys.argv[1]

all_files = []
for file in glob.glob(parent_dir + "/*scatterPlot.pdf"):
	file = file.split("/")[-1]
	all_files.append(file[:-20])

all_files = sorted(set(all_files))

#print(all_files)
#exit(0)


f = open(parent_dir + "/all_SNV_plots.tex" , 'w')

f.write("\\documentclass[15pt]{article}\n \\usepackage[english]{babel}\n \\usepackage[utf8x]{inputenc}\n \\usepackage{graphicx}\n \\usepackage[margin=1in]{geometry}\n \\usepackage[font=Large]{caption}\n \\begin{document}\n\n")

for file in all_files:

	patient = file.split("_")[-2]
	sample = file.split("_")[-1]

	f.write("\\centering\\section*{Patient = " + "{}".format(patient) + " ; sample = " + "{}".format(sample) + " ; primer = M1}\n \\vspace*{0.5in}\n\\centering\n \\begin{figure}[h]\n \\centering\n \\textbf{\\Large{Repeat 1 \\hspace*{1.9in} Repeat 2}}\\par\\medskip\n \\includegraphics[width=0.45\\textwidth]{./" + file + "_M1A_scatterPlot.pdf}\n \\includegraphics[width=0.45\\textwidth]{./" + file + "_M1B_scatterPlot.pdf}\n \\end{figure}\n\n \\includegraphics[width=0.8\\textwidth]{./" + file + "_M1_replicate_frequencies.png}\\\\ \n \\clearpage \n\n")
	f.write("\\centering\\section*{Patient = " + "{}".format(patient) + " ; sample = " + "{}".format(sample) + " ; primer = M2}\n \\vspace*{0.5in}\n\\centering\n \\begin{figure}[h]\n \\centering\n \\textbf{\\Large{Repeat 1 \\hspace*{1.9in} Repeat 2}}\\par\\medskip\n \\includegraphics[width=0.45\\textwidth]{./" + file + "_M2A_scatterPlot.pdf}\n \\includegraphics[width=0.45\\textwidth]{./" + file + "_M2B_scatterPlot.pdf}\n \\end{figure}\n\n \\includegraphics[width=0.8\\textwidth]{./" + file + "_M2_replicate_frequencies.png}\\\\ \n \\clearpage \n\n")



f.write("\\end{document}\n")
f.close()