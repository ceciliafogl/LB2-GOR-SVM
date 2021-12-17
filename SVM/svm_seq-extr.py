import sys
import os

if __name__=='__main__':
	input_file = sys.argv[1]		#blind set ids
	svm_pred_file = sys.argv[2]		#prediction per each model
	with open(input_file) as file_in:
		svm_prediction = open(svm_pred_file).read().splitlines()
		index = 0
		length = 0
		for ids in file_in:
			ids = ids.rstrip()
			out = open(ids+".seqs", "w+")
			O = '>'
			O += ids + "\n"
			z = ''
			#fasta_file = open('/Users/ceciliafoglini/Desktop/Lab2/project/fasta/'+ids+'.fasta').read().splitlines()
			fasta_file = open('/Users/ceciliafoglini/Desktop/Lab2/project/150Random_DSSPs_fasta/'+ids+'.dssp.fasta').read().splitlines() #for blind set
			seq = fasta_file[1]
			length += len(seq)
			for i in seq:
				while index < length:
					#print(index, length)
					if svm_prediction[index] == '1':
						O += 'H'
						index += 1
					elif svm_prediction[index] == '2':
						O += 'E'
						index += 1
					elif svm_prediction[index] == '3':
						O += '-'
						index += 1
					else:
						print('index error', index)
			out.write(O)
			out.close()