
import numpy as np 
import sys
import argparse


# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('--input_one', help='', type=str)
parser.add_argument('--input_two', help='', type=str)
parser.add_argument('--input_three', help='', type=str)
parser.add_argument('--output', help='', type=str)
args = parser.parse_args()


# Read in data files
position , p_val , raw_freq = np.loadtxt(args.input_one , unpack=True , usecols=(2,5,6))
ref_base , var_base = [] , []

for line in open(args.input_one , 'r').readlines():
	fields = line.replace('   ',' ').replace('  ' , ' ' ).split(" ")
	ref_base.append(fields[3])
	var_base.append(fields[4])
	#ref_base.append(line.split('   ')[1])
	#var_base.append(line.split('   ')[2].split(' ')[0])

#ref_base , var_base = np.genfromtxt(args.input_one , unpack=True , usecols=(3,4))

n_tst_fw , cov_tst_fw , n_tst_bw , cov_tst_bw , n_ctrl_fw = np.loadtxt(args.input_two , unpack=True , skiprows=1 , usecols=(2,3,4,5,6))
#cov_tst_fw , n_tst_bw , cov_tst_bw , n_ctrl_fw , cov_ctrl_fw , n_ctrl_bw , cov_ctrl_bw = np.loadtxt(args.input_two , unpack=True , usecols=(1,2,3,4,5,6,7))
cov_ctrl_fw , n_ctrl_bw , cov_ctrl_bw = np.loadtxt(args.input_three , unpack=True , skiprows=1 , usecols=(1,2,3))


# Compute "shifted" variant frequencies
#shifted_var_freq = (( n_tst_fw + n_tst_bw )/( cov_tst_fw + cov_tst_bw )) - (( n_ctrl_fw + n_ctrl_bw )/( cov_ctrl_fw + cov_ctrl_bw ))
#print(n_tst_fw)
#print(n_tst_bw)
#print(cov_tst_fw)
#print(cov_tst_bw)
shifted_var_freq = (( n_tst_fw + n_tst_bw )/( cov_tst_fw + cov_tst_bw ))


filtered_out = []

# Filtering 
f = open(args.output + '.METRICS.dat' , 'w')

if not isinstance(n_tst_bw, np.float64):
	for i in range(len(n_tst_bw)):

		# Filter on raw number of calls foreach variant
		if ((n_tst_bw[i] + n_tst_fw[i]) <= 10):
			filtered_out.append(i)
			f.write("Removed entry {}{}{} due to small number of raw calls\n".format(int(position[i]) , ref_base[i] , var_base[i]))

		# Filter on shifted_var_freq value, which accounts for coverage at particular position
		#if (shifted_var_freq[i] < 0.0025):
		#	filtered_out.append(i)
		#	f.write("Removed entry {} due to frequency below threshold\n".format(i))

elif isinstance(n_tst_bw, np.float64):
	
	# Filter on raw number of calls foreach variant
	if ((n_tst_bw + n_tst_fw) <= 10):
		filtered_out.append(0)
		f.write("Removed entry {}{}{} due to small number of raw calls\n".format(int(position[i]) , ref_base[i] , var_base[i]))

f.close()




# Open output file 
if not isinstance(n_tst_bw, np.float64):

	g = open(args.output , 'w')

	for i in range(len(shifted_var_freq)):
		if (i in filtered_out):
			continue
		else:
			g.write("{}{}{} {:1.10f} {}\n".format(int(position[i]) , ref_base[i] , var_base[i] , shifted_var_freq[i] , p_val[i]))

	g.close()



	g = open(args.output + ".noGermline.dat" , 'w')
	for i in range(len(shifted_var_freq)):

		# If VAF in control is greater than 1%, count mutation as germline and do not write to file
		if (i in filtered_out):
			continue

		if (( n_ctrl_fw[i] + n_ctrl_bw[i] )/( cov_ctrl_fw[i] + cov_ctrl_bw[i] ) <= 0.01):
		 	g.write("{}{}{} {:1.10f} {}\n".format(int(position[i]) , ref_base[i] , var_base[i] , shifted_var_freq[i] , p_val[i]))


	g.close()



else:

	g = open(args.output , 'w')

	if (len(filtered_out) == 0):
		g.write("{}{}{} {:1.10f} {}\n".format(int(position) , ref_base , var_base , shifted_var_freq , p_val))


	g.close()


	g = open(args.output + ".noGermline.dat" , 'w')

	if (len(filtered_out) == 0) and (( n_ctrl_fw + n_ctrl_bw )/( cov_ctrl_fw + cov_ctrl_bw ) <= 0.01):
		g.write("{}{}{} {:1.10f} {}\n".format(int(position) , ref_base , var_base , shifted_var_freq , p_val))

	g.close()















