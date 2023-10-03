#!/usr/bin/env python

# type: ignore
# pylint: skip-file

#################################################
# Transforms oligomap output to SAM format and	
# keeps only best alignments.			
#						
# input:					
# - the output of Oligomap	
# - nh filter value		
# output:					
# - stdout: a SAM file with the best		
#   alignments for each read *(e.g. all	3	
#   alignments with 0 errors, or all alignments	
#   with 1 error for each read
# 
#  Paula Iborra. Zavolan Lab.
#  Adated version of Alessandro Crippa script.
#################################################


import sys
import re
from argparse import ArgumentParser, RawTextHelpFormatter


### ARGUMENTS ### PI, modified June 2019. 

parser = ArgumentParser(
    description="Oligomap output to SAM. NH filter applicable."
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '-n','--nhfilter',
    help="Add NH tag to output, remove reads that contain more aligments than given NH value (with min error).",
    type=int, 
    )
parser.add_argument(
    '-i', '--input',
    help="Input File. The output of oligomap mapping.",
    required=True, 
    )

args = parser.parse_args()


readSeqs = {}		# dictionary of reads
filt = {}		# dict with all filtered reads (heavy multimappers NH > 100) by error. {'seqName':error} (error with which it has been discarted)
seqToMinError = {}
nh = {}			# nh dictionary per read evaluated


def addReadToList(d,nh,minerr,seqName,flag,target,positionInTarget,errors,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString):
	
	if len(d) == 0: # if empty dictionary
		if seqName not in list(minerr.keys()): # means: seqName NOT found previously
			d[seqName] = []
			minerr[seqName] = errors
			nh[seqName] = 1
			d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])
		else: # means: seqName found previously and filtered by NH ---- empty dicty due to d.clear() 
			if errors < minerr[seqName]: # Aligment found with lower error than the ones stored in minerr dict
				d[seqName] = []
				minerr[seqName] = errors
				nh[seqName] = 1
				d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])

	else: # if dictionary not empty
		if seqName == list(d.keys())[0]:	#same seqName as the one stored in dict found
			# check the errors fo this new seqName to include it or not
			if errors == minerr[seqName]:	# if seqName errors is equal to the one stored in minerr -> keep 
				nh[seqName] += 1 # increase nh +1 
				if args.nhfilter:	# if NH filter 
					if nh[seqName] <= (args.nhfilter):
						d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])
					else:	# if after adding +1 to nh, the total is > than the filter value. seqName discarded. 
						# Clear d dictonary
						d.clear()
						sys.stderr.write("Filtered by NH | Read %s | Errors = %s \n"%(seqName,errors))
						
				else: 	# no NH filtering, keep all reads including heavy multimappers
					d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])

			elif errors < minerr[seqName]: 	# error minor to the one previously stored.
					sys.stderr.write("Filtered by ERROR | Read %s | Errors = %s \n"%(seqName,minerr[seqName])) 
					minerr[seqName] = min(errors, minerr[seqName])	# Update minor error.
					d[seqName] = [] 	# Create empy list for seqName (removed aligmants previously stored with higher error).
					nh[seqName] = 1 	# NH starts at 1
					d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])
			else: 	# if seqName errors is > than the error stored, discard 
				pass
		
		elif seqName != list(d.keys())[0]: #new seqName found
			# PRINT to the output file the previous seqName stored in dict
			for i in d.keys():
				sys.stderr.write("Printed read %s | Errors = %s | NH = %s \n"%(i, minerr[i], nh[i]))
				for al in d[i]:
					nhtag= ('NH:i:'+str(nh[i]))
					al.append(nhtag)
					print('\t'.join([str(x) for x in al]))

			# Clean dictonaries before saving the current seqName being evaluated
			d.clear()
			nh.clear()
			minerr.clear()

			# Restart the dictionaries with the new seqName
			d[seqName] = []
			minerr[seqName] = errors
			nh[seqName] = 1
			d[seqName].append([seqName,flag,target,positionInTarget,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString])
	


def readInput(fi):
	i = 0
	#read oligomap results
	while True:
		line1 = fi.readline()	#e.g. seq_61 (21 nc) 21..1  chr1  1..21
		if not line1: break	#[check if file is empty]
		line2 = fi.readline()	#e.g. chr1
		if not line2: break	#[pointless check]
		line3 = fi.readline()   #e.g. errors: 1 orientation: -
		if not line3: break	#[pointless check]
		line4 = fi.readline()   #e.g. CAAGCAGAAGACGGCATACGC	->input sequence
		if not line4: break	#[pointless check]
		line5 = fi.readline()   #e.g. |||||||||||||||||||| 
		if not line5: break	#[pointless check]
		line6 = fi.readline()   #e.g. CAAGCAGAAGACGGCATACGA	->target sequence
		if not line6: break	#[pointless check]
		line7 = fi.readline()   #e.g. ""
		i += 1
		#we got a result
		if line1[0].isdigit():	#out readNames start with a number
			#parse lines
			line1 = re.sub( '\s+', ' ', line1).strip().split()
			line3 = re.sub( '\s+', ' ', line3).strip().split()
			
			seqName = line1[0].strip()
			target = line2.strip()
			positionInTarget = line1[5].split('.')[0].strip()
			#position2InTarget= line1[5].split('.')[2].strip()
			errors = line3[1].strip()
			strand = line3[3].strip()
			sequence = line4.strip()
			referenceSeq = line6.strip()
			mapq='255'
			rnext = '*'
			pnext= '0'
			tlen='0'
			qual='*'	

			#flag
			if strand == "+":
				flag = "0"
			else:
				flag = "16"
			#define cigar string
			if errors == '0':
			#perfect match, cigar string is just the number of nucleotides + "M"
				cigarStr = str(len(sequence)) + "M"
				matchingString = str(len(sequence))
				editDistance = "NM:i:0"
			else:
			#cigar and mismatch strings are built here
				if line5[0] == " ":					#if the first nucleotide is mutated
					cigarStr= "1M"
				elif line5[-2] == " ":					#if the last nucleotide is mutated
					cigarStr= str(len(line5)-2) +"M"
				else:
					cigarStr= str(line5.strip().index(' ')) + "M"	#if any other nucleotide is mutated
				editDistance = "NM:i:1"
				#depending on where in the read a deletion,insertion or mutation occurs, we create ad hoc mismatch strings.
				#case 1: deletion in the seq
				if '-' in sequence:
					indelerr = "1D"
					if line5[0] == " ":		#the 1st nt is deleted
						cigarStr = indelerr + str(len(sequence)-1) + "M"
						matchingString = "^" + referenceSeq[0] + str(len(sequence)-1)
					elif line5[-2] == " ":	#the last nt is deleted ([-2] because [-1] is "\n")
						cigarStr = str(len(sequence)-1) + "M" + indelerr
						matchingString = str(len(sequence)-1) + "^" + referenceSeq[len(sequence)-1] + "0"	#"0" is required if the last entry is ^[something]
					else:			#deletion occurs "in" the read
						tmp = cigarStr
						cigarStr = cigarStr + indelerr + str(line5.strip().count('|')-int(tmp[:-1])) + "M"
						matchingString = str(line5.strip().index(' ')) + "^" + referenceSeq[int(tmp[:-1])] + str(len(sequence) - line5.strip().index(' ') -1)
				#case 2: insertion in the seq
				elif '-' in referenceSeq:
					indelerr = "1I"
					matchingString = str(len(sequence))
					if line5[0] == " ":		#the 1st nt is inserted
						cigarStr = indelerr + str(len(sequence)-1) + "M"
					elif line5[-2] == " ":	#the last nt is inserted
						cigarStr = str(len(sequence)-1) + "M" + indelerr
					else:			#addition occurs "in" the read
						tmp = cigarStr
						cigarStr = cigarStr + indelerr + str(line5.strip().count('|')-int(tmp[:-1])) + "M"
				#case 3: single point mutation
				else:
					cigarStr = str(len(sequence)) + "M"
					if line5[0] == " ":		#the 1st nt is inserted
						matchingString = referenceSeq[0] + str(len(sequence)-1)
					elif line5[-2] == " ":	#the last nt is inserted
						matchingString = str(len(sequence)-1) + referenceSeq[len(sequence)-1]
					else:			#mutation occurs "in" the read
						matchingString = str(line5.strip().index(' ')) + referenceSeq[line5.strip().index(' ')] + str(len(sequence) - line5.strip().index(' ') -1)
					#indelerr = referenceSeq[int(cigarStr)] #this info is not used in the cigar string
				sequence = re.sub( '-', '', sequence).strip()

			matchingString = ('MD:Z:'+matchingString)		
			#addReadToList(readSeqs, readFilt, readNh, seqName,flag,target,positionInTarget,errors,cigarStr,sequence,editDistance,matchingString)
			sys.stderr.write("Record: %i | Sequence: %s \n"%(i,seqName))
			addReadToList(readSeqs, nh, seqToMinError, seqName,flag,target,positionInTarget,errors,mapq,cigarStr,rnext,pnext,tlen,sequence,qual,editDistance,matchingString)
			

fi = open(args.input, "r")	#read file
sys.stderr.write("###########################\nSTART READING...\n###########################\n")
readInput(fi)			#process
fi.close()

#print last aligments stored in dict
if len(readSeqs) != 0: 
	for i in readSeqs.keys():
		sys.stderr.write("Printed read %s | Errors = %s | NH = %s \n"%(i, seqToMinError[i], nh[i]))
		for al in readSeqs[i]:
			nhtag= ('NH:i:'+str(nh[i]))
			al.append(nhtag)
			print('\t'.join([str(x) for x in al]))

sys.stderr.write("SUCCESSFULLY FINISHED.")
sys.exit()
