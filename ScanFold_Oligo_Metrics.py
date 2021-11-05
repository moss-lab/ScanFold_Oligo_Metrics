#ScanFold_Oligo_Metrics

"""
By Warren B. Rouse, Moss Lab, Iowa State University

Usage:python ScanFold_Oligo_Metrics.py  FinalPartners.txt  Sequence.fasta  ScanFold.log  Outfile  IntegerWindowSize SequenceID

#Script requires the ScanFold Final Partners file, input fasta file, ScanFold log file, an outfile name without an extension, oligo size as an integer value, and a sequence ID to find the average of the average z-scores, mfe, and ed across a user defined window size; number of pairs in the window; number of base pairs per nucleotide in the ScanFold window; reverse complement of the sequence in the window; the range of nucleotides; and the sequence ID to be analyzed

"""

import os
import sys
import numpy as np

zs =[]							#List of avgzscores
zavg = []						#List of 19nt average of avgzscores					
mfe = []						#List of avgMFEs
avgmfe =[]						#List of 19nt average of avgMFEs
ed = []							#List of avgEDs	
avged = []						#List of 19nt average of avgEDs
bpi = []						#List of bpi
bpj = []						#List of bpj
paired = []						#List of bpi and bpj that are paired as 1 and unpaired as 0
pairedness = []					        #List of concatenated number of pairs in each 19nt window
seq = []						#List of individual reverse complemented nts
seqlist = []					        #List of all sliced reversed complemented nucleotides
revcomp = []						#List of reversed reverse complemented nts
new_seq = []						#List of reversed reverse complement sequence strings to list new_seq
nums = []						#List of BPs per nt
BPs = []						#List of sliced BP per nt
BPlist =[]						#List of sliced BP per nt without fringe cases

partners = sys.argv[1]					#Defining Final partners file needed for to run the program as partners
sequence = sys.argv[2]					#Defining fasta file needed for to run the program as sequence
log = sys.argv[3]					#Defining log file needed for to run the program as log
output_name = sys.argv[4]				#Defining user input file needed for to run the program as output_name
window = sys.argv[5]				    	#Defining user input for window size being analyzed (i.e. 19, 20, 21, etc.)
ID = sys.argv[6]					#Defining user input for transcript ID 

with open(partners , 'r') as final:			#Opens Final Partners file as final
	lines = final.readlines()[1:]			#Read all lines in the input file except the first header line
	for line in lines:				#For loop to read values from final partners file
		data = line.split()			#Lines into separate data columns
		bpi.append((int(data[1])))		#Appends bpi values to list bpi
		bpj.append((int(data[2])))		#Appends bpj values to list bpj
		zs.append((float(data[4])))		#Appends avgzscore values to list zs
		mfe.append((float(data[3])))		#Appends mfe values to list mfe
		ed.append((float(data[5])))		#Appends ed values to list ed
		if data[1]==data[2]:			#If statement for bpi values equal to bpj values
			paired.append(0)		#Appends value of 0 to list paired when bpi and bpj are equal
		else:					#Else statement for bpi values not equal to bpj values
			paired.append(1)		#Appends value of 1 to list paired when bpi and bpj are not equal
		#print(paired)

i = 0								#Defining position i
j = (int(window))						#Defining position j
for x in paired:						#For loop to iterate through bpi and bpj values in paired list
	summedpair = (sum(paired[i:j]))				#Create a list of summed values from position i to j in window size of interest
	if len(paired[i:j])==int(window):			#If statement to define when to write data to list pairedness
		pairedness.append((summedpair))			#Appending full length (19nt) data to list pairedness
	else:							#Else statement to define when not to write data to list pairedness
		pass						#Do not write to list if data is not full length (19nt)
	i += 1							#Defining position i + 1 to move down list by 1
	j += 1							#Defining position j + 1 to move down list by 1
#print(pairedness)

i = 0								#Defining position i
j = (int(window))						#Defining position j
for x in zs:							#For loop to iterate through zavg values in list
	avglist = round(np.mean(zs[i:j]),2)			#Create a list of mean values from position i to j and round to 3 decimal places
	if len(zs[i:j])==int(window):				#If statement to define when to write data to list zavg
		zavg.append((avglist))				#Appending full length (19nt) data to list zavg
	else:							#Else statement to define when not to write data to list zavg
		pass						#Do not write to list if data is not full length (19nt)
	i += 1							#Defining position i + 1 to move down list by 1
	j += 1							#Defining position j + 1 to move down list by 1
#print(zavg)
											
i = 0								#Defining position i
j = (int(window))						#Defining position j
for x in mfe:							#For loop to iterate through avgmfe values in list
	mfelist = round(np.mean(mfe[i:j]),2)			#Create a list of mean values from position i to j and round to 3 decimal places
	if len(mfe[i:j])==int(window):				#If statement to define when to write data to list avgmfe
		avgmfe.append((mfelist))			#Appending full length (19nt) data to list avgmfe
	else:							#Else statement to define when not to write data to list avgmfe
		pass						#Do not write to list if data is not full length (19nt)
	i += 1							#Defining position i + 1 to move down list by 1
	j += 1							#Defining position j + 1 to move down list by 1
#print(len(avgmfe))

i = 0								#Defining position i
j = (int(window))						#Defining position j
for x in ed:							#For loop to iterate through avgED values in list
	edlist = round(np.mean(ed[i:j]),2)			#Create a list of mean values from position i to j and round to 3 decimal places
	if len(ed[i:j])==int(window):				#If statement to define when to write data to list avged
		avged.append((edlist))				#Appending full length (19nt) data to list avged
	else:							#Else statement to define when not to write data to list avged
		pass						#Do not write to list if data is not full length (19nt)
	i += 1							#Defining position i + 1 to move down list by 1
	j += 1							#Defining position j + 1 to move down list by 1
#print(len(avged))

def reverse_complement(dna):									#Defining the function reverse_complement
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 'A': 'U'}			#Defining what the complement nucleotides are
    return ''.join([complement[base] for base in dna[::-1]])					#Join individual complemented nucleotides in a string and reverse the list

with open(sequence , 'r') as final:				#Opens Final Partners file as final
	lines = final.readlines()[1:]				#Read all lines in the input file except the first header line
	for line in lines:					#For loop to read values from fasta file
		data = reverse_complement(line)			#Creating a reverse complement of the sequence
		revcomp.append(data[::-1])			#Append reversed reverse complement to list revcomp
	#print(revcomp)
	for x in revcomp[0]:					#For loop to iterate through the reverse complement list	
		seq.append(x)						#Append individual positions to the seq list
	#print(seq)

i = 0							#Defining position i
j = (int(window))					#Defining position j
for x in seq:						#For loop to iterate through position in reverse complement list seq
	allseqs = seq[i:j]				#Create a list of sliced sequences
	i += 1						#Defining position i + 1 to move down list by 1
	j += 1						#Defining position j + 1 to move down list by 1
	if len(allseqs)==int(window):			#If statement to define when to write data to list seqlist
		seqlist.append((allseqs))		#Appending full length (19nt) data to list seqlist
	else:						#Else statement to define when not to write data to list seqlist
		pass					#Do not write to list if data is not full length (19nt)
#print(seqlist)
								
for x in seqlist:					#For loop to iterate through sequences in list
	seqstring = ''					#Setting seqstring as a list
	for y in x:					#For loop to iterate through positions in the sequence
		seqstring += y 				#Add each posistion to string
	new_seq.append(seqstring)			#Append sequence strings to list new_seq
#print(new_seq)

with open(log , 'r') as final:				#Opens ScanFold log file as final
	lines = final.readlines()			#Read all lines in the input file
	for line in lines:				#For loop to read values from log file
		data = line.split()			#Lines into separate data columns
		if 'BPs' in line:			#If statement to pull out only lines with string BPs in it
			#print(line)
			nums.append(data[10])		#Appending data BPs per nt to the list nums
	#print(nums)

i = 0								#Defining position i
j = (int(window))						#Defining position j
for x in nums:							#For loop to iterate through BP per nt in list nums
	BPs = nums[i:j]						#Defining BPs as a list of sliced BP per nt
	i += 1							#Defining position i + 1 to move down list by 1
	j += 1							#Defining position j + 1 to move down list by 1
	if len(BPs)==int(window):				#If statement to define when to write data to list BPs
		BPlist.append((BPs))				#Appending full length (19nt) data to list BPs
	else:							#Else statement to define when not to write data to list BPs
		pass   						#Do not write to list if data is not full length (19nt)
#print(BPlist)

with open(f"{output_name}.txt", "w") as file: 																																																																		#Open out file to write to
	i = 1																																																																											#Defining position i for writing to outfile
	j = int(window)																																																																									#Defining position j for writing to outfile
	file.write(f'average z-score, MFE, ED; #ofBPs; #ofBP/nt; reverse complement of oligo; range of nucleotides being analyzed in a tiled {window} nt window; and Sequence ID\n')																																	#Writing description of file header
	file.write('zavg\tavgMFE\tavgED\t#ofBPs\tpos1BP/nt\tpos2BP/nt\tpos3BP/nt\tpos4BP/nt\tpos5BP/nt\tpos6BP/nt\tpos7BP/nt\tpos8BP/nt\tpos9BP/nt\tpos10BP/nt\tpos11BP/nt\tpos12BP/nt\tpos13BP/nt\tpos14BP/nt\tpos15BP/nt\tpos16BP/nt\tpos17BP/nt\tpos18BP/nt\tReverse Complement\tntRange\tSequence ID\n')			#Writing header and going to next line in out file 
	for z in range(0,len(zavg)):																																																																					#For loop to write values to out file
		if z == (len(zavg) - 1):																																																																					#If statement for writing last value to file
			file.write(f'{zavg[z]}\t{avgmfe[z]}\t{avged[z]}\t{pairedness[z]}\t{"	".join(map(str, BPlist[z]))}\t{new_seq[z]}\t{i}-{j}\t{ID}')																																										#Write last value and corresponding nt number without adding extra empty line
		else:  																																																																										#Else statement for writing all values other than the last
			file.write(f'{zavg[z]}\t{avgmfe[z]}\t{avged[z]}\t{pairedness[z]}\t{"	".join(map(str, BPlist[z]))}\t{new_seq[z]}\t{i}-{j}\t{ID}\n')																																									#Write all values other than the last with corresponding nt number and add new line
			i += 1																																																																									#Defining position i + 1 to write the nt range correspoing to the avg value
			j += 1																																																																									#Defining position j + 1 to write the nt range correspoing to the avg value
file.close()																																																																										#Close the file after writing it
