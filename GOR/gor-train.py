import sys
import os
import numpy as np

#Extract seq profile

def matrix_pssm(filename):
	seq_prof = []
	for j in filename:
		L = j.split()
		try: L[0]
		except: pass
		else:
			if L[0].isdigit():
				seq_prof.append(L[22:-2])
	x = np.array(seq_prof, dtype=np.float64)
	x /= 100
	if x.sum() != float(0):
		return(x), True
	else:
		return(None), False

def profile_matrices(profile1, dssp1, H1, E1, C1, TOT1, SS1):
	pad = np.zeros((8,20))
	padding = np.vstack((pad,profile1,pad))
	Lp=len(padding)
	dssp_seq = dssp1.read().splitlines()[1]
	dssp_ok='-'*8 + dssp_seq + '-'*8
	for e in range(8,Lp-8):
		j=e-8
		k=e+8
		if dssp_ok[e] == 'H':
			H1 += padding[j:k+1]
			SS1[0] += 1
		if dssp_ok[e] == 'E':
			E1 += padding[j:k+1]
			SS1[1] += 1
		if dssp_ok[e] == '-':
			C1 += padding[j:k+1]
			SS1[2] += 1
		TOT1 += padding[j:k+1]
	return(H1,E1,C1,TOT1,SS1)


if __name__ == '__main__':
	fileid = sys.argv[1]
	H = np.zeros(shape=(17,20))
	E = np.zeros(shape=(17,20))
	C = np.zeros(shape=(17,20))
	TOT = np.zeros(shape=(17,20))
	SS = np.zeros(3)

	with open(fileid) as filein:
		for id in filein:
			id=id.rstrip()
			profile_file = '/Users/ceciliafoglini/Desktop/Lab2/project/jpred_pssm/' + id + '.fasta.pssm'
			dssp_file = '/Users/ceciliafoglini/Desktop/Lab2/project/dssp/' + id + '.dssp'
			try:
				prof=open(profile_file)
				dssp=open(dssp_file)
			except: continue
			else:
				profile, ret = matrix_pssm(prof)
				if ret:
					H, E, C, TOT, SS = profile_matrices(profile, dssp, H, E, C, TOT, SS)	#SS= alfa,beta,coil
	N_res=SS[0]+SS[1]+SS[2]
	H /= N_res
	E /= N_res
	C /= N_res
	TOT /= N_res
	SS /= N_res

	
	np.save('train-out_H',H)
	np.save('train-out_E',E)
	np.save('train-out_C',C)
	np.save('train-out_TOT',TOT)
	np.save('train-out_SS',SS)
