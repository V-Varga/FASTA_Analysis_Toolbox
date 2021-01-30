#!/bin/python
"""
Title: fasta_an.py
Date: 2020-10-28
Authors: Virág Varga

Description:
	This program takes a FASTA file and performs various analyses on the file:
		- Quality-Checking: This function is integrated into the other functions. It ensures
			that the input FASTA file is of proper quality for analysis.
		- 6-Frame Translation: This function translates all sequences in the FASTA file
			into all 6 possible reading frames. It then evaluates which protein sequence
			is most likely to be the real reading frame, based on the lengths of the generated
			protein sequences.
		- Nucleotide Frequency Calculation: This function calculates the frequency of the 4
			standard nucleotide bases present in each sequence in the FASTA file. The GC
			content of the file is also calculated. A pdf file containing pie charts displaying
			graphical representations of these frequencies is created.
		- BLAST @ NCBI: This function uses Biopython to BLAST the FASTA file against the
			NCBI database. Results are outputted to a specified file. It is recommended that
			this file be in .xml format.

List of functions:
	qualityCheck
	readFrame1
	readFrame2
	readFrame3
	readFrame4
	readFrame5
	readFrame6
	likelyProtein
	translateFasta
	GCcontent
	FASTAblast

List of standard and non-standard modules used:
	argparse
	os
	matplotlib (requires installation)
	BioPython (requires installation)

Procedure:
	1. Assignment of command-line arguments; initial quality-checking of input FASTA file.
	2. Function designated by flag is executed:
		a. -t: translate
			Executed when called.
			Uses the standard codon table to translate the nucleotide sequence into all 6 possible
			reading frames for the protein product.
			Highlights the longest protein sequence(s) as the most likely actual protein
			product of the nucleotide sequence.
		b. -c: GCcontent
			Executed when called.
			Uses counters to calculate the frequency of the 4 standard nucleotide bases in
			each sequence. These frequencies are used to calculate the GC content and
			produce a pdf containing pie charts showing the relative frequencies of the
			nucleotide bases for each sequence.
		c. '-b': FASTAblast
			Executes when called.
			Performs a BLASTn on the input FASTA sequence using the NCBI database. Produces
			an xml-formatted file of the results.
	3. Outputs from functions are written out to the designated files.

Known bugs and limitations:
	- Files for which the id/header lines do not start with a '>' character cannot be used as inputs.
	- No non-standard nucleotide characters can be accepted by the program for translation.
	- If illegal characters are included in the file name/path, the code will not run, and will return
		an error.
	- Only nucleotide FASTA files can be used as inputs; amino acid FASTA are not accepted.
	- Currently this program is only able to handle 1 flag at a time. The program must be run multiple
		times with the unique flags and specified outfiles to produce the different outputs.
	- Owing to the way the GCcontent function works (flag: '-c'), the pdf output of the image results is
		named as specified in the function, not by the user designating the name of the .pdf file in the
		command line call. To change the pdf file name, the user must alter that portion of the script.
	- Both the input and output files must have specified names in order for the program to run.
	- This program will not work properly for cases where the sequence associated with a FASTA header is
		not in one line, but instead split over multiple lines.

Usage:
	./fasta_an.py [-h] [-t] [-c] [-b] [-v] input_file output_file
	OR
	python fasta_an.py [-h] [-t] [-c] [-b] [-v] input_file output_file

This script was written for Python 3.8.5, in Spyder 4.
"""

#################################   ARGPARSE   #######################################
import argparse
#the argparse module allows for a single program script to be able to carry out a variety of specified functions
#this can be done with the specification of unique flags for each command


parser = argparse.ArgumentParser(description = 'This program can perform various analyses on a standard DNA nucleotide FASTA file input.')
#The most general description of what this program can do is defined here


parser.add_argument(
	'-t', '--translate',
	action='store_true',
	help = 'This argument will translate the nucleotide FASTA into the protein sequences resulting from the 6 possible reading frames, and evaluate the most likely protein product(s).'
	)
	#the '-t' flag will call the translateFasta() function to be executed on the input FASTA file
	#the results will be written out to the specified output file
parser.add_argument(
	'-c', '--gc_content',
	action='store_true',
	help = 'This argument will calculate the frequencies of the 4 standard DNA nucleotide bases for each sequence in the FASTA file, along with the GC content.'
	)
	#the '-c' flag will call the GCcontent() function to be executed on the input FASTA file
	#the text-based results will be written out to the specified output file
	#the graphical results will be written into a pdf file whose name is specified in the script (ie. not user-defined)
parser.add_argument(
	'-b', '--BLAST',
	action='store_true',
	help = 'This argument will BLAST the input FASTA file against the NCBI database.'
	)
	#the '-b' flag will call the FASTAblast() function to be executed on the input FASTA file
	#the results will be written to the specified output file in XML format
	#if the user specifies their output file to be in .xml format, they will be able to open the results directly in a browser
parser.add_argument(
	#'-i', '--input',
	#the above line of code is left in as further clarification of this argument
	dest='input_file',
	metavar='INPUT_FILE',
	type=argparse.FileType('r')
	)
	#this portion of code specifies that the program requires an input file, and it should be opened for reading ('r')
parser.add_argument(
	#'-o', '--output',
	#the above line of code is left in as further clarification of this argument
	dest='output_file',
	metavar='OUTPUT_FILE',
	type=argparse.FileType('w')
	)
	#this portion of code specifies that the program requires an output file, and it should be opened for writing ('w')
parser.add_argument(
	'-v', '--version',
	action='version',
	version='%(prog)s 1.0'
	)
	#This portion of the code specifies the version of the program; currently 1.0
	#The user can call this flag ('-v') without specifying input and output files


args = parser.parse_args()
#this command allows the program to execute the arguments in the flags specified above



###############################   QUALITY CHECKING   #################################

#The quality-checking function used throughout this program is defined here
#This function is called inside of all other flag-defined functions in this program to ensure FASTA quality before analysis


def qualityCheck(infile):
	#quality-checking will only be performed on the input file
	first_line = infile.readline().strip()
	#the contents of the first line, without the '\n' character, are assigned to variable 'first_line'
	if not first_line.startswith('>'):
		raise Exception('Input file not in proper FASTA format.')
		#if the first line of the FASTA file does not begin with a '>' character the code will raise this exception
	for line in infile:
		#the fucntion will read through the FASTA file line by line
		acceptChars = set(['A', 'T', 'G', 'C'])
		#a set of acceptable characters is defined here
		#non-standard DNA nucleotide bases are not accepted for analysis by this program
		if not line.startswith('>'):
			#only lines that do not start with a '>' character (ie. non-header lines) will be checked for nucleotide base contents
			seq = line.strip().upper()
			#the sequence is defined in variable 'seq' as containing the FASTA sequence without the '\n' end-line character
			#the bases are also turned uppercase, to avoid issues with case-matching
			seqChars = set(seq)
			#the members of the sequence defined in variable 'seq' are converted into a set defined in variable 'seqChars'
			for chars in seqChars:
				#this portion of the code will iterate through the characters in set 'seqChars'
				if chars not in acceptChars:
					raise Exception('Non-standard nucleotides included in FASTA file.')
					#if the sequence contains characters that are not a standard DNA nucleotide base, this exception will be raised
				else:
					return infile
					#the input file will be returned for use in the rest of the function



#################################   TRANSLATE   ######################################

#This portion of the code deals with the translation of the nucleotide FASTA sequences into their 6 possible reading frames
#flag '-t' can be used to call it from the command line


table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
#The codon translation table used throughout translation is defined in the dictionary 'table'


#The functions below translate the 6 possible reading frames

#Translating the 3 forward-direction reading frames:

def readFrame1(seq):
	#calculation of the 1st reading frame
	#uses the nucleotide sequence as input
	protein1 =""
	#variable 'protein1' is defined as an empty string
	for i in range(0, len(seq), 3):
		#the code will begin translation starting at the 1st nucleotide base in the sequence
		#translation will continue for the length of the sequence, in sets of 3 nucleotides (ie. codons)
		if (len(seq)-(i+2)) > 2:
			#translation will continue as long as the length of the sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide sequence portion
			#and lasting through the 3rd nucleotide in the sequence portion
			if table[codon] == '_':
				protein1 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the sequence will then stop
			else:
				protein1 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein1
	#the string 'protein1' will be returned
	#this contains the 1st translated protein sequence


def readFrame2(seq):
	#calculation of the 2nd reading frame
	#uses the nucleotide sequence as input
	protein2 =""
	#variable 'protein2' is defined as an empty string
	for i in range(1, len(seq), 3):
		#the code will begin translation starting at the 2nd nucleotide base in the sequence
		#translation will continue for the length of the sequence, in sets of 3 nucleotides (ie. codons)
		if (len(seq)-(i+2)) > 2:
			#translation will continue as long as the length of the sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide sequence portion
			#and lasting through the 3rd nucleotide in the sequence portion
			if table[codon] == '_':
				protein2 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the sequence will then stop
			else:
				protein2 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein2
	#the string 'protein2' will be returned
	#this contains the 2nd translated protein sequence


def readFrame3(seq):
	#calculation of the 3rd reading frame
	#uses the nucleotide sequence as input
	protein3 =""
	#variable 'protein3' is defined as an empty string
	for i in range(2, len(seq), 3):
		#the code will begin translation starting at the 3rd nucleotide base in the sequence
		#translation will continue for the length of the sequence, in sets of 3 nucleotides (ie. codons)
		if (len(seq)-(i+2)) > 2:
			#translation will continue as long as the length of the sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide sequence portion
			#and lasting through the 3rd nucleotide in the sequence portion
			if table[codon] == '_':
				protein3 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the sequence will then stop
			else:
				protein3 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein3
	#the string 'protein3' will be returned
	#this contains the 3rd translated protein sequence


#Translating the 3 reverse-direction reading frames:

def readFrame4(rev_seq):
	#calculation of the 4th reading frame (1st reverse)
	#uses the reversed nucleotide sequence as input
	protein4 =""
	#variable 'protein4' defined as an empty string
	for i in range(0, len(rev_seq), 3):
		#the code will begin translation starting at the 1st nucleotide base in the reversed sequence
		#translation will continue for the length of the reversed sequence, in sets of 3 nucleotides (ie. codons)
		if (len(rev_seq)-(i+2)) > 2:
			#translation will continue as long as the length of the reversed sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = rev_seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide reversed sequence portion
			#and lasting through the 3rd nucleotide in the reversed sequence portion
			if table[codon] == '_':
				protein4 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the reversed sequence will then stop
			else:
				protein4 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein4
	#the string 'protein4' will be returned
	#this contains the 4th translated protein sequence


def readFrame5(rev_seq):
	#calculation of the 5th reading frame (2nd reverse)
	#uses reversed nucleotide sequence as input
	protein5 =""
	#variable 'protein5' defined as an empty string
	for i in range(1, len(rev_seq), 3):
		#the code will begin translation starting at the 1st nucleotide base in the reversed sequence
		#translation will continue for the length of the reversed sequence, in sets of 3 nucleotides (ie. codons)
		if (len(rev_seq)-(i+2)) > 2:
			#translation will continue as long as the length of the reversed sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = rev_seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide reversed sequence portion
			#and lasting through the 3rd nucleotide in the reversed sequence portion
			if table[codon] == '_':
				protein5 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the reversed sequence will then stop
			else:
				protein5 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein5
	#the string 'protein5' will be returned
	#this contains the 5th translated protein sequence


def readFrame6(rev_seq):
	#calculation of the 6th reading frame (3rd reverse)
	#uses reversed nucleotide sequence as input
	protein6 =""
	#variable 'protein6' defined as an empty string
	for i in range(2, len(rev_seq), 3):
		#the code will begin translation starting at the 1st nucleotide base in the reversed sequence
		#translation will continue for the length of the reversed sequence, in sets of 3 nucleotides (ie. codons)
		if (len(rev_seq)-(i+2)) > 2:
			#translation will continue as long as the length of the reversed sequence is >2 nucleotides
			#(ie. can still be translated into a 3-nucleotide codon)
			codon = rev_seq[i:i + 3]
			#variable 'codon' defines the translation codon as starting from the 1st nucleotide in the 3-nucleotide reversed sequence portion
			#and lasting through the 3rd nucleotide in the reversed sequence portion
			if table[codon] == '_':
				protein6 += table[codon]
				break
				#when a codon translates to a '_' character (ie. a stop codon), the '_' will be printed into the protein sequence
				#translation for this reading frame of the reversed sequence will then stop
			else:
				protein6 += table[codon]
				#if the codon specified is not a stop codon ('_' character), then the amino acid will be printed into the string
				#translation will continue
		else:
			break
			#once there are fewer than 3 nucleotides left, translation will stop
	return protein6
	#the string 'protein6' will be returned
	#this contains the 6th translated protein sequence


#Comparing the translated protein sequence lengths
def likelyProtein(frameList):
	#uses a list containing the 6 translated protein sequence strings as input
	max_length = max(frameList, key=len)
	#finds the longest sequence(s) in list 'frameList'
	max_length = len(max_length)
	#calculates the length of the longest sequence(s) in list 'frameList'
	max_prot = [i for i in frameList if len(i) == max_length]
	#creates a list 'max_prot' containing all sequences in list 'frameList' that have the length of the longest sequence
	return max_prot
	#returns list 'max_prot'


#Using all of the functions defined above to do the overall translation for the entire FASTA file:

def translateFasta(infile, outfile):
	#uses the input file as the source of data, and the output file to write the results into
	qualityCheck(infile)
	#checks the quality of the input FASTA file before analysis is performed
	for line in infile:
		#reads through the FASTA file line by line
		if line.startswith('>'):
			#identifies the header lines
			header = line.lstrip('>')
			#variable 'header' contains the header line without the '>' character
			outfile.write(header)
			#writes the header line into the output file
		else:
			#identifies the sequence line that follows the header line
			seq = line.replace('\n', '').upper()
			#variable 'seq' contains the sequence associated with the header
			#the end-line character '\n' is removed and all bases are converted to uppercase to prevent issues with case-matching
			prot1 = readFrame1(seq)
			#the product of function readframe1() is saved in variable 'prot1'
			prot2 = readFrame2(seq)
			#the product of function readframe2() is saved in variable 'prot2'
			prot3 = readFrame3(seq)
			#the product of function readframe3() is saved in variable 'prot3'
			outfile.write('Forward Protein Sequences:' + '\n' + prot1 + '\n' + prot2 + '\n' + prot3 + '\n')
			#the translations of the 3 forward-direction reading frames are written out to the output file
			rev_seq = seq[::-1]
			#variable 'rev_seq' contains the reverse of the nucleotide sequence from variable 'seq'
			prot4 = readFrame4(rev_seq)
			#the product of function readframe4() is saved in variable 'prot4'
			prot5 = readFrame5(rev_seq)
			#the product of function readframe5() is saved in variable 'prot5'
			prot6 = readFrame6(rev_seq)
			#the product of function readframe6() is saved in variable 'prot6'
			outfile.write('Reverse Protein Sequences:' + '\n' + prot4 + '\n' + prot5 + '\n' + prot6 + '\n')
			#the translations of the 3 reverse-direction reading frames are written out to the output file
			frameList = [prot1, prot2, prot3, prot4, prot5, prot6]
			#list 'frameList' is constructed, containing the 6 translated protein sequences
			likelyseq = likelyProtein(frameList)
			#variable 'likelyseq' contains the list of sequences with the greatest length
			likelyseq = map(lambda x:x+'\n', likelyseq)
			#the contents of variable 'likelyseq' are converted into a format that can be written out to a file
			outfile.write('Most Likely Protein Sequence(s):' + '\n')
			outfile.writelines(likelyseq)
			outfile.write('\n')
			#the most likely protein sequence(s) are written out to the output file

'''
It should be noted that the way the protein sequences are printed, '_' will be printed if there is a stop codon.
As a result, it is conceivable that 2 protein sequences would have the same "length" if one ends with a '_' and another
has a translation that goes until the end of the nucleotide sequence.
The user should be able to determine for themselves if this is the case, by examing the output file.
'''

if args.translate:
	translateFasta(args.input_file, args.output_file)
	#when the '-t' flag is used in the command-line call, this command is executed
	#the contents of the input nucleotide FASTA file be translated into the 6 possible reading frames
	#the output file of user-determined name will be structured as defined above



###########################   GC content   #############################################

#This portion of the code deals with calculating the frequencies of the nucleotide bases in each sequence in the input FASTA file
#flag '-c' can be used to call it from the command line

import matplotlib.pyplot as plt
#this module allows the constuction of the pie chart products
from matplotlib.backends.backend_pdf import PdfPages
#this module allows the pie chart products to be printed to a pdf file


def GCcontent(infile, outfile):
	#the contents of the input FASTA file will be used as the source of data, and the analysis results will be written out to the output file
	qualityCheck(infile)
	#checks the quality of the input FASTA file before analysis is performed
	img_count = 0
	#counter 'img_count' will be used to make sure that each pie chart of the nucleotide frequencies is saved and printed to the result pdf file
	count_dict = {}
	#dictionary 'count_dict' will contain the frequencies of the nucleotide bases
	header_dict = {}
	#dictionary 'header_dict' will contain the headers associated with the nucleotide sequences and sequence counters
	for line in infile:
		#reads through the FASTA file line by line
		if line.startswith('>'):
			#identifies the header lines by the '>' character
			header = line.lstrip('>')
			#variable 'header' contains the header line without the '>' character
			outfile.write(header)
			#the header line is written out into the output file
		else:
			#identifies the sequence lines
			seq = line.replace('\n', '').upper()
			#variable 'seq' contains the sequence line without the end-line character '\n'
			#nucleotide bases are all converted to capital letters to avoid issues with case-matching
			G_count = 0
			#counter 'G_count' will count the instances of 'G's in the sequence; it is reset to 0 for each new sequence
			C_count = 0
			#counter 'C_count' will count the instances of 'C's in the sequence; it is reset to 0 for each new sequence
			A_count = 0
			#counter 'A_count' will count the instances of 'A's in the sequence; it is reset to 0 for each new sequence
			T_count = 0
			#counter 'T_count' will count the instances of 'T's in the sequence; it is reset to 0 for each new sequence
			for i in seq:
				#reads through the sequence base by base
				if i == 'G':
					G_count += 1
					#if the nucleotide base is a 'G' +1 is added to the counter 'G_count'
				if i == 'C':
					C_count += 1
					#if the nucleotide base is a 'C' +1 is added to the counter 'C_count'
				if i == 'A':
					A_count += 1
					#if the nucleotide base is a 'A' +1 is added to the counter 'A_count'
				if i == 'T':
					T_count += 1
					#if the nucleotide base is a 'T' +1 is added to the counter 'T_count'
			seq_len = len(seq)
			#variable 'seq_len' contains the total calculated length of the sequence
			GC_cont = (G_count + C_count) / seq_len * 100
			#variable 'GC_cont' contains the calculated GC content of the sequence
			outfile.write('The nucleotide percentages are:' + '\n')
			outfile.write('G: ' + str(G_count) + '/' + str(seq_len) + '\t' + 'C: ' + str(C_count) + '/' + str(seq_len) + '\n')
			outfile.write('A: ' + str(A_count) + '/' + str(seq_len) + '\t' + 'T: ' + str(T_count) + '/' + str(seq_len) + '\n')
			#the frequencies of the nucleotide bases are written out to the output file
			outfile.write('The GC content of the sequence is: ' + '\n')
			outfile.write("{:.3f}".format(GC_cont) + '%' + '\n\n')
			#the GC content of the sequence is written out to the output file
			count_dict[img_count] = [G_count, C_count, A_count, T_count]
			#dictionary 'count_dict' is filled with the contents of counter 'img_count' as the key
			#and a list containing the contents of the 4 nucleotide counters as the associated value
			header_dict[img_count] = header
			#dictionary 'header_dict' is filled with the header for the sequence associated with the counter 'img_count' associated with that sequence
			img_count += 1
			#counter 'img_count' is increased by +1 for the next iteration of the loop
	pp = PdfPages('base_freqs.pdf')
	#variable 'pp' will contain the opened pdf file that will contain the pie charts showing the nucleotide frequencies for each sequence
	for img_count, GCAT_values in count_dict.items():
		#iterates through the dictionary 'count_dict'
		plt.figure(img_count)
		#a figure (pie chart) will be created for each individual 'img_count' key
		labels = 'G', 'C', 'A', 'T'
		#the labels for the "slices" in the pie chart are defined here
		plt.pie(GCAT_values, labels=labels, autopct='%1.1f%%', startangle=90)
		#the contents of the pie chart are defined as the members of the list contained in the value associated with the defined 'img_count' key
		#% (percent) values are displayed in the chart for the frequency of the nucleotide bases
		plt.title(header_dict[img_count])
		#each pie chart is given the title of the header associated with that sequence
		pp.savefig()
		#the figure is saved to the pdf output file
	#plt.show()
	#the above line of code is left in so that if the user wishes to run this code in Spyder, or a similar visualizer
	#the code will then print the results to the visualizer's window
	pp.close()
	#closes the pdf result file

'''
As noted in the header portion of this program, the name of the output pdf file containing the image results is pre-defined in the above code.
If the user wishes to change the name of the output pdf file, they must do so by altering that portion of the script.
Note that this limits the usability of this function on multiple FASTA files contained in the same directory.
'''

if args.gc_content:
	GCcontent(args.input_file, args.output_file)
	#when the '-c' flag is used in the command-line call, this command is executed
	#the contents of the input nucleotide FASTA file will be used to calculate the nucleotide frequencies for each individual sequence
	#the output file of user-determined name will be structured as defined above
	#the output file containing the image results will be named as defined in the code above



##############################   BLAST @ NCBI  ###########################################################################

#This portion of the code will BLAST the contents of the input FASTA file against the NCBI database, and return a results file in XML format
#flag '-b' can be used to call it from the command line

from Bio.Blast import NCBIWWW
#this module uses BioPython to conduct the BLAST


def FASTAblast(input_file, output_file):
	#the input and output files should be defined by the user in the command line
	#it is recommended that the output file be in .xml format, to allow it to be opened directly in a browser
	qualityCheck(input_file)
	#checks the quality of the input FASTA file before analysis is performed
	fastafile = input_file
	#defines the FASTA file to be BLAST-ed as the input file
	result_handle = NCBIWWW.qblast('blastn', 'nt', fastafile.read())
	#performs a BLASTn of the input nucleotide FASTA file against the NCBI database
	output_file.write(result_handle.read())
	#writes the results to the user-defined output file in XML format


if args.BLAST:
	FASTAblast(args.input_file, args.output_file)
	#when the '-b' flag is used in the command-line call, this command is executed
	#the contents of the input nucleotide FASTA file will be BLAST-ed against the NCBI database
	#the output file of user-determined name will be structured as defined above
