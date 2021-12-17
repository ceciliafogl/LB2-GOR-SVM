## Generate input files for SVM ##

import sys
import os
import numpy as np
np.set_printoptions(threshold=np.inf)

#extract sequence profile

def matrix_pssm(filename, number, ID):
	n = number
	seq_prof = []
	for j in filename:
		L = j.split()
		try: L[0]
		except: pass
		else:
			if L[0].isdigit():
				seq_prof.append(L[22:-2])
	n += len(seq_prof)
	x = np.array(seq_prof, dtype=np.float64)
	x /= 100
	if np.sum(x) != float(0):
		return(x), True , n
	else:
		return(None), False, n

def SVM_input(profile1, dssp1, n):
	P = len(profile1)
	pad = np.zeros((8,20))
	SVM_profile = np.zeros(shape=(P,340))
	padding = np.vstack((pad,profile1,pad))
	Lp = len(padding)
	dssp_seq = dssp1.read().splitlines()[1]
	dssp_ok='-'*8 + dssp_seq + '-'*8
	for e,i in zip(range(8,Lp-8), range(0,P)):
		j=e-8
		k=e+8
		window = padding[j:k+1]
		SVM_line = window.flatten()
		SVM_profile[i] += SVM_line
	class_column = np.zeros(shape=(P,1))
	SVM_profile_w_class = np.append(class_column,SVM_profile,axis=1)
	SVM_profile_list_w_class = SVM_profile_w_class.tolist()
	for e in range(0,P):						
		if dssp_ok[e] == 'H':
			SVM_profile_list_w_class[e][0] = 1
		if dssp_ok[e] == 'E':
			SVM_profile_list_w_class[e][0] = 2
		if dssp_ok[e] == '-':
			SVM_profile_list_w_class[e][0] = 3
	for line in SVM_profile_list_w_class:
		for indice1 in range(0,341):
			line[indice1] = indice1,':', line[indice1]			
		filtered_line = [i for i in line if i[2] > 0]				
		filtered_line[0] = filtered_line[0][2]
		print(filtered_line)

if __name__ == '__main__':
	fileid = sys.argv[1]
	no_pssm = open("blind_no_pssm.txt", "w+")
	NP = ""
	N = 0
	no_hits = open("blind_no_hits", "w+")
	NH = ""
	with open(fileid) as filein:
		for ids in filein:
			ids=ids.rstrip()
#			profile_file = '/Users/ceciliafoglini/Desktop/Lab2/project/jpred_pssm/' + ids + '.fasta.pssm'		
#			dssp_file = '/Users/ceciliafoglini/Desktop/Lab2/project/dssp/' + ids + '.dssp'
			profile_file = '/Users/ceciliafoglini/Desktop/Lab2/project/Blind-Set_pssm/'+ids+'.dssp.fasta.pssm'	#blind prediction
			dssp_file = '/Users/ceciliafoglini/Desktop/Lab2/project/150Random_DSSPs_dssp/'+ids+'.ss.dssp'
			try:
				prof=open(profile_file)
				dssp=open(dssp_file)
			except:
				NP += ids+'\n'
				continue
			else:
				profile, ret, N = matrix_pssm(prof, N, ids)
				if ret:
					SVM_input(profile, dssp, N)

				else:
					print(ids+'\n')
	no_pssm.write(NP)
	no_pssm.close()
	no_hits.write(NH)
	no_hits.close()
