
import numpy as np 
import sys
import argparse


# Parse command line arguments
parser = argparse.ArgumentParser()
#parser.add_argument('-I', help='Input file with raw coverage data')
parser.add_argument('--input_one', help='', type=str)
parser.add_argument('--input_two', help='', type=str)
parser.add_argument('--output', help='', type=str)
args = parser.parse_args()


# Read in data files
position , p_val , raw_freq , n_tst_fw = np.loadtxt(args.input_one , unpack=True , usecols=(2,5,6,8))
ref_base , var_base = [] , []

for line in open(args.input_one , 'r').readlines():
	fields = line.replace('   ',' ').replace('  ' , ' ' ).split(" ")
	ref_base.append(fields[3])
	var_base.append(fields[4])
	#ref_base.append(line.split('   ')[1])
	#var_base.append(line.split('   ')[2].split(' ')[0])

#ref_base , var_base = np.genfromtxt(args.input_one , unpack=True , usecols=(3,4))

cov_tst_fw , n_tst_bw , cov_tst_bw , n_ctrl_fw , cov_ctrl_fw , n_ctrl_bw , cov_ctrl_bw = np.loadtxt(args.input_two , unpack=True , skiprows=1 , usecols=(1,2,3,4,5,6,7))


# Compute "shifted" variant frequencies
shifted_var_freq = (( n_tst_fw + n_tst_bw )/( cov_tst_fw + cov_tst_bw )) - (( n_ctrl_fw + n_ctrl_bw )/( cov_ctrl_fw + cov_ctrl_bw ))



# Open output file 
f = open(args.output , 'w')

for i in range(len(shifted_var_freq)):
	f.write("{}{}{} {:1.10f} {}\n".format(int(position[i]) , ref_base[i] , var_base[i] , shifted_var_freq[i] , p_val[i]))

f.close()