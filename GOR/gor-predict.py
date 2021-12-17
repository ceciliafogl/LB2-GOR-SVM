import sys
import os
import numpy as np
import math


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
		return x, True
	else:
		return None, False


def information(helix, beta, coil, total, secstr):  #information function
	Mh=np.array(total*secstr[0])	#denominator
	Me=np.array(total*secstr[1])
	Mc=np.array(total*secstr[2])
	Mh1=helix/Mh			#fraction
	Me1=beta/Me
	Mc1=coil/Mc
	Ih_r=np.log(Mh1)		#matrix
	Ie_r=np.log(Me1)
	Ic_r=np.log(Mc1)
	return(Ih_r, Ie_r, Ic_r)



def prediction(profile1, Ihr1, Ier1, Icr1):
	pad = np.zeros((8,20))
	padding = np.vstack((pad,profile1,pad))
	Lp=len(padding)
	dssp=''
	ss = ['H', 'E', '-']
	for e in range(8,Lp-8):
		j=e-8
		k=e+8
		W=np.array(padding[j:k+1])	
		prediction=[]
		Ih=np.sum(W*Ih_r)
		Ie=np.sum(W*Ie_r)
		Ic=np.sum(W*Ic_r)
		prediction=[Ih,Ie,Ic]
		SS_pred=max(prediction)
		dssp += ss[prediction.index(SS_pred)]
	return dssp


if __name__ == '__main__':
	HELIX=np.load('train-out_H.npy')
	BETA=np.load('train-out_E.npy')
	COIL=np.load('train-out_C.npy')
	TOTAL=np.load('train-out_TOT.npy')
	SECSTR=np.load('train-out_SS.npy')
	blind_fileid=sys.argv[1]

	with open(blind_fileid) as filein:
		for id in filein:
			id=id.rstrip()
			with open('SS-prediction_'+id+'.dssp', 'w+') as output:
				profile_file = '/Users/ceciliafoglini/Desktop/Lab2/project/Blind-Set_pssm/' + id + '.dssp.fasta.pssm'
			#	profile_file = '/Users/ceciliafoglini/Desktop/Lab2/project/jpred_pssm/' + id + '.fasta.pssm' #CV
				try:
					prof=open(profile_file, "r")
				except: continue
				else:
					profile, ret  = matrix_pssm(prof)
					if ret:
						Ih_r, Ie_r, Ic_r = information(HELIX, BETA, COIL, TOTAL, SECSTR)
						dssp = prediction(profile, Ih_r, Ie_r, Ic_r) #I(S;R)=out di information(h,e,c,tot,secstr)
					output.write('>'+id+'\n'+dssp+'\n')
