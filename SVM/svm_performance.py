import sys
import os
import numpy as np
import math


def confusion_matrix(matrice_totale,seq1, seq2):
	residues_seq1 = len(seq1)
	residues_seq2 = len(seq2)
	conf_mat = [[0,0,0],[0,0,0],[0,0,0]]

	for res in (range(residues_seq1)):
		if seq1[res] == seq2[res]:
			if seq1[res] == 'H':
				conf_mat[0][0] += 1
			elif seq1[res] == 'E':
				conf_mat[1][1] += 1
			elif seq1[res] == '-':
				conf_mat[2][2] += 1
		else:
			if seq1[res] == 'H':
				if seq2[res] == 'E':
					conf_mat[1][0] += 1
				elif seq2[res] == '-':
	  				conf_mat[2][0] += 1
			elif seq1[res] == 'E':
				if seq2[res] == 'H':
					conf_mat[0][1] += 1
				elif seq2[res] == '-':
					conf_mat[2][1] += 1
			elif seq1[res] == '-':
				if seq2[res] == 'H':
					conf_mat[0][2] += 1
				elif seq2[res] == 'E':
					conf_mat[1][2] += 1

	npConf_Mat=np.array(conf_mat)
	return(npConf_Mat)


def performance(c_mat):

	tot_c_mat = (c_mat[0][0]+c_mat[1][1]+c_mat[2][2]+c_mat[0][1]+c_mat[0][2]+c_mat[1][0]+c_mat[1][2]+c_mat[2][0]+c_mat[2][1]) #is it a np array
	ACC = (c_mat[0][0]+c_mat[1][1]+c_mat[2][2])/tot_c_mat

# with respect to H
	Ch = c_mat[0][0]
	Nh = c_mat[1][1] + c_mat[1][2] + c_mat[2][1] + c_mat[2][2]
	Oh = c_mat[0][1] + c_mat[0][2]
	Uh = c_mat[1][0] + c_mat[2][0]
	try:
		SEN_h = (Ch / (Ch+Uh))
		PPV_h = (Ch / (Ch+Oh))
	except ZeroDivisionError:
		SEN_h = 0
		PPV_h = 0

	MCC_h = ((Ch*Nh)-(Oh*Uh)) / math.sqrt((Ch+Oh)*(Ch+Uh)*(Nh+Oh)*(Nh+Uh))
# with respect to E
	Ce = c_mat[1][1]
	Ne = c_mat[0][0] + c_mat[0][2] + c_mat[2][0] + c_mat[2][2]
	Oe = c_mat[1][0] + c_mat[1][2]
	Ue = c_mat[0][1] + c_mat[2][1]
	try:
		SEN_e = Ce / (Ce+Ue)
		PPV_e = Ce / (Ce+Oe)
	except ZeroDivisionError:
		SEN_e = 0
		PPV_e = 0
	MCC_e = ((Ce*Ne)-(Oe*Ue)) / math.sqrt((Ce+Oe)*(Ce+Ue)*(Ne+Oe)*(Ne+Ue))
# with respect to C
	Cc = c_mat[2][2]
	Nc = c_mat[0][0] + c_mat[0][1] + c_mat[1][0] + c_mat[1][1]
	Oc = c_mat[2][0] + c_mat[2][1]
	Uc = c_mat[0][2] + c_mat[1][2]
	try:
		SEN_c = Cc / (Cc+Uc)
		PPV_c = Cc / (Cc+Oc)
	except ZeroDivisionError:
		SEN_c = 0
		PPV_c = 0

	MCC_c= ((Cc*Nc)-(Oc*Uc)) / math.sqrt((Cc+Oc)*(Cc+Uc)*(Nc+Oc)*(Nc+Uc))
####################

	print('Q3 = '+str(ACC))
	print('SEN_h = '+str(SEN_h), 'PPV_h = '+str(PPV_h),'MCC_h = '+str(MCC_h))
	print('SEN_e = '+str(SEN_e), 'PPV_e = '+str(PPV_e),'MCC_e = '+str(MCC_e))
	print('SEN_c = '+str(SEN_c), 'PPV_c = '+str(PPV_c),'MCC_c = '+str(MCC_c))





if __name__ == '__main__':
	print('CONFUSION MATRIX\n')
	IdsPredFile = sys.argv[1]
#	CV_set = sys.argv[2]
	conf_mat_TOTALE = np.zeros((3, 3))
#	print(conf_mat_TOTALE)
	with open(IdsPredFile) as filein:
		for id in filein:
			id = id.rstrip()
#			print('>'+id)
			observed_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/150Random_DSSPs_dssp/'+id+'.ss.dssp'
			predicted_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/SVM/for_project/final_predicted_sequences/'+id+'.seqs'
			#predicted_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/SVM/final_predicted_blind_sequences/'+id+'.seqs'
			with open(observed_SecStr) as observed_SS, open(predicted_SecStr) as predicted_SS:
				for line1, line2 in zip(observed_SS, predicted_SS):
					if len(line1) != len(line2):
						print(id)
						print(len(line1))
						print(len(line2))
						continue
					if line1[0] == '>': continue
					if line2[0] == '>': continue	#AND oppure OR condition
					else:
						conf_mat_TOTALE += confusion_matrix(conf_mat_TOTALE, line1, line2)
#						print(len(line1))
#						print(len(line2))
#						print("\n")
#			break
	print(conf_mat_TOTALE)
	PERFORM_MAT = performance(conf_mat_TOTALE)
