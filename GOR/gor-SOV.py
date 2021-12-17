## Segment OVerlap Index for the prediction ##

import sys
import os
import numpy as np
import math
import re
from statistics import mean


def segment_finder(seq_o, seq_p):
	seqs = [seq_o, seq_p]
	frag = {'H':[[],[]], 'E':[[], []], '-':[[], []]}
	rr = ['H', 'E', '-']

	for c in rr:
		my_regex = re.escape(c) + r"+"
		for seq, num in zip(seqs,range(2)):
			idx = 0
			while c in seq[idx:]:
				sr = seq[idx:]
				sing_frag = []
				x = re.search(my_regex, sr)

				for i in range(x.span()[0]+idx, x.span()[1]+idx):
					sing_frag.append(i)

				frag[c][num].append(sorted(sing_frag))
				idx += x.span()[1]

	return frag['H'][0], frag['H'][1], frag['E'][0], frag['E'][1], frag['-'][0], frag['-'][1]


def SOV(seq_o, seq_p, SegFind_out):
	frag = {'H':[SegFind_out[0],SegFind_out[1]], 'E':[SegFind_out[2],SegFind_out[3]], 'C':[SegFind_out[4],SegFind_out[5]]}
	minov_residues_H = []
	minov_H = []
	minov_residues_E = []
	minov_E = []
	minov_residues_C = []
	minov_C = []

	maxov_residues_H = []
	maxov_H = []
	maxov_residues_E = []
	maxov_E = []
	maxov_residues_C = []
	maxov_C = []

	N_H = 0
	N_E = 0
	N_C = 0

	length_o_frags = []
	length_p_frags = []
	for ss,frags in frag.items():								
		obs_frags = frags[0]
		pred_frags = frags[1]

		for obs_frag in obs_frags:
			n_res_in_obs_frag = len(obs_frag)
			if ss == 'H':
				N_H += n_res_in_obs_frag
			elif ss == 'E':
				N_E += n_res_in_obs_frag
			elif ss == 'C':
				N_C += n_res_in_obs_frag

			for pred_frag in pred_frags:
				intersection = set(obs_frag) & set(pred_frag)
				INT = sorted(intersection)
				union = set(obs_frag) | set(pred_frag)
				UN = sorted(union)
				if len(INT) > 0:
					if ss == 'H':
						minov_H.append(len(INT))
						minov_residues_H.append(INT)
						maxov_H.append(len(UN))
						maxov_residues_H.append(UN)
					if ss == 'E':
						minov_E.append(len(INT))
						minov_residues_E.append(INT)
						maxov_E.append(len(UN))
						maxov_residues_E.append(UN)
					if ss == 'C':
						minov_C.append(len(INT))
						minov_residues_C.append(INT)
						maxov_C.append(len(UN))
						maxov_residues_C.append(UN)


	len_S_o_H = []
	len_S_p_H = []
	half_len_S_o_H = []
	half_len_S_p_H = []

	len_S_o_E = []
	len_S_p_E = []
	half_len_S_o_E = []
	half_len_S_p_E = []

	len_S_o_C = []
	len_S_p_C = []
	half_len_S_o_C = []
	half_len_S_p_C = []

	for ss,frags in frag.items():
		obs_frags = frags[0]
		pred_frags = frags[1]
		for obs_frag in obs_frags:
			for pred_frag in pred_frags:
				intersection = set(obs_frag) & set(pred_frag)
				INT = sorted(intersection)
				if len(INT) > 0:
					if ss == 'H':
						len_S_p_H.append(len(pred_frag))
						half_len_S_p_H.append(float((len(pred_frag))/2))
						len_S_o_H.append(len(obs_frag))
						half_len_S_o_H.append(float((len(obs_frag))/2))
					if ss == 'E':
						len_S_p_E.append(len(pred_frag))
						half_len_S_p_E.append(float((len(pred_frag))/2))
						len_S_o_E.append(len(obs_frag))
						half_len_S_o_E.append(float((len(obs_frag))/2))
					if ss == 'C':
						len_S_p_C.append(len(pred_frag))
						half_len_S_p_C.append(float((len(pred_frag))/2))
						len_S_o_C.append(len(obs_frag))
						half_len_S_o_C.append(float((len(obs_frag))/2))



	N_residues = [N_H, N_E, N_C]						
	l_maxov = [maxov_H, maxov_E, maxov_C]					
	l_minov = [minov_H, minov_E, minov_C]					
	l_len_o = [len_S_o_H, len_S_o_E, len_S_o_C]				
	l_len_p = [len_S_p_H, len_S_p_E, len_S_p_C]				
	l_half_len_o = [half_len_S_o_H, half_len_S_o_E, half_len_S_o_C]		
	l_half_len_p = [half_len_S_p_H, half_len_S_p_E, half_len_S_p_C]		

	dlt_H = []
	dlt_E = []
	dlt_C =[]
	deltas = [dlt_H, dlt_E, dlt_C]
	rr = ['H', 'E', 'C']

	SOV_H_frags = []
	SOV_E_frags = []
	SOV_C_frags = []

	SOV_H = []
	SOV_E = []
	SOV_C = []


	for i in range(len(N_residues)):
		N = N_residues[i]
		for k in range(len(l_maxov[i])):
			delta_arg = [l_maxov[i][k] - l_minov[i][k], l_minov[i][k], l_half_len_o[i][k], l_half_len_p[i][k]]
			delta = min(delta_arg)
			if rr[i] == 'H':
				dlt_H.append(delta)
			if rr[i] == 'E':
				dlt_E.append(delta)
			if rr[i] == 'C':
				dlt_C.append(delta)

			SOV_frag = ((l_minov[i][k]+delta)/l_maxov[i][k])*l_len_o[i][k]

			if rr[i] == 'H':
				SOV_H_frags.append(SOV_frag)
			if rr[i] == 'E':
                                SOV_E_frags.append(SOV_frag)
			if rr[i] == 'C':
                                SOV_C_frags.append(SOV_frag)


	return(SOV_H_frags, SOV_E_frags, SOV_C_frags, N_residues)

def SOV_calculation(sov_h_frags, sov_e_frags, sov_c_frags, N_ss):
	sov_h_prot = 100*(1/N_ss[0])*(sum(sov_h_frags))
	sov_e_prot = 100*(1/N_ss[1])*(sum(sov_e_frags))
	sov_c_prot = 100*(1/N_ss[2])*(sum(sov_c_frags))

	return(sov_h_prot, sov_e_prot,sov_c_prot)


if __name__ == '__main__':
	print('SOV:\n')
	tot_prot = 0
	IdsPredFile = sys.argv[1]
	with open(IdsPredFile) as filein:
		SOV_H_sum = 0
		SOV_E_sum = 0
		SOV_C_sum = 0
		for id in filein:
			tot_prot += 1
			id = id.rstrip()
			observed_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/150Random_DSSPs_dssp/' + id + '.ss.dssp'
			#observed_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/dssp/' + id + '.dssp' #cv
			predicted_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/GOR/gor_blind_predictions/SS-prediction_' + id + '.dssp'
			#predicted_SecStr = '/Users/ceciliafoglini/Desktop/Lab2/project/GOR_cv/test4/pred_seq/SS-prediction_' + id + '.dssp' #
			with open(observed_SecStr,'r') as observed_SS:
				with open(predicted_SecStr, 'r') as predicted_SS:
					observed_SS = observed_SS.read().splitlines()
					predicted_SS = predicted_SS.read().splitlines()
					H_obs, H_pred, E_obs, E_pred, C_obs, C_pred = segment_finder(observed_SS[1], predicted_SS[1])
					SF_out = [H_obs, H_pred, E_obs, E_pred, C_obs, C_pred]
					SOV_H_frags, SOV_E_frags, SOV_C_frags, N_secstr = SOV(observed_SS[1], predicted_SS[1], SF_out)
					try:
						SOV_H_prot, SOV_E_prot, SOV_C_prot = SOV_calculation(SOV_H_frags, SOV_E_frags, SOV_C_frags, N_secstr)
					except:
						continue
					SOV_H_sum += SOV_H_prot
					SOV_E_sum += SOV_E_prot
					SOV_C_sum += SOV_C_prot

		SOV_H_tot = SOV_H_sum/tot_prot
		SOV_E_tot = SOV_E_sum/tot_prot
		SOV_C_tot = SOV_C_sum/tot_prot

		print(SOV_H_tot, SOV_E_tot, SOV_C_tot)
