#!/usr/bin/python3

import os
import sys
from Bio.Seq import Seq


######################## RETRIEVE SEQ STARTS ###########################

def retrieve(nt_i, nt_f, genome_dir, genome_revcomp):
	'''
	It retrieves a subsequence from the sequence or from the reverse complementary sequence.
	It compares the first and the last nt of the subsequence to determine from which chain it should retrieve that subseq.
	Returns the subsequence.
	'''
	len_genome = len(genome_dir)
	if nt_f > nt_i:
		return genome_dir[nt_i - 1:nt_f]
	else:
		return genome_revcomp[len_genome - nt_i:len_genome - nt_f + 1]

def retrieve_run(genome_fasta_fname, positions_list_fname, rp_option = False, rp_from = False, rp_to = False):
	'''
	Checks if the user files exist and launches the main function.
	It alse checks for the optional parameters
	'''
	genome_fname_ok = try_open_file(genome_fasta_fname)
	positions_fname_ok = try_open_file(positions_list_fname)
	if rp_option:
		if rp_to <= rp_from:
			print('Relative position of first nucletide has to be smaller than last nucleotide')
			quit()

	retrieve_main(genome_fname_ok, positions_fname_ok, rp_option, rp_from, rp_to)

def relativeposition(orf_nt_i, orf_nt_f, rp_option, rp_from, rp_to):
	'''
	OPTIONAL FUNCTION - RELATIVE SEQUENCE RETRIEVE
	It calculates the first and last nt positions of the sequence in relation to a relative gene nt.
	It returns the first and the last nt position.
	(first_nt, last_nt)
	'''
	if orf_nt_f > orf_nt_i:
		orf_orientation = 'dir'
	else:
		orf_orientation = 'rev'
	
	if rp_option == '-start':
		nt_i = orf_nt_i
		nt_f = orf_nt_i
	elif rp_option == '-stop':
		nt_i = orf_nt_f
		nt_f = orf_nt_f
	elif rp_option == '-startstop':
		nt_i = orf_nt_i
		nt_f = orf_nt_f
	
	if orf_orientation == 'dir':
		return nt_i + rp_from, nt_f + rp_to
	elif orf_orientation == 'rev':
		return nt_i - rp_from, nt_f - rp_to

def retrieve_main(genome_fasta_fname, positions_list_fname, rp_option, rp_from, rp_to):
	'''
	Opens the txt file with the positions of every subsequence and the genome FASTA file.
	Loops through it's lines.
	Saves the subsequences in a FASTA file.
	'''
	genome_fasta = open(genome_fasta_fname)
	positions_list = open(positions_list_fname)
	if rp_option:
		file_out = open(genome_fasta_fname.replace('.FASTA', f'{rp_option}_{str(rp_from)}_{str(rp_to)}.FASTA'), 'a')
	if not rp_option:
		file_out = open(genome_fasta_fname.replace('.FASTA', '-sequences.FASTA'), 'a')

	for line in genome_fasta:
		line = line.strip()
		if line.startswith('>'):
			continue
		genome_dir = Seq(line)
		break #if its a multiFASTA, it only reads the first sequence
	genome_revcomp = genome_dir.reverse_complement()
	genome_len = len(genome_dir)

	for i, line in enumerate(positions_list):
		line = line.strip().split()
		if not line:
			continue
		try: #for missing data points or wrong data type
			gene_name, nt_i, nt_f = line[0], int(line[1]), int(line[2])

			if rp_option:
				nt_i, nt_f = relativeposition(nt_i, nt_f, rp_option, rp_from, rp_to)

			if nt_i < 0 or nt_f < 0:
				print(f'[Error] - Negative position.')
				print(f'  Sequence {gene_name} in file {positions_list_fname} skipped.')
				print(f'  Line {i + 1}: {gene_name} {nt_i} {nt_f}')
				continue
			if nt_i > genome_len or nt_f > genome_len:
				print(f'[Error] - Position bigger that genome size.')
				print(f'  Sequence {gene_name} in file {positions_list_fname} skipped.')
				print(f'  Line {i + 1}: {gene_name} {nt_i} {nt_f}')
				continue
		except: #index or value error
			print(f'[Error] - Wrong data point.')
			print(f'  Sequence {line[0]} in file {positions_list_fname} skipped.')
			print(f'  Line {i + 1}: {(" ").join(line)}')
			continue

		file_out.write(f'>{gene_name}\n{retrieve(nt_i, nt_f, genome_dir, genome_revcomp)}\n')

	file_out.close()
	genome_fasta.close()
	positions_list.close()


######################## RETRIEVE SEQ ENDS ###########################

def try_open_file(file_name):
	'''
	Tries to open the users file. 
	If it fails, gives an error message and quits.
	Returns file_name.
	'''
	try:
		file_ready = open(file_name)
		file_ready.close()
	except:
		print(file_name, 'not found')
		quit()
	return file_name

def all_files(file_extension = '.FASTA'):
	'''
	Makes a list with all the file names in that folder ending with file_extension.
	By default it looks for .FASTA files.
	'''
	files = [x_file for x_file in os.listdir(os.curdir) if x_file.endswith(file_extension)]
	if not files: #no files in that folder with the file_extension extension
		print('No', file_extension, 'files found')
		quit()
	files.sort()
	return files


def print_help():
	help_list = ['Retrieve the sequence of every ORF from a genome or retrieve sequences relative to the first, last or both nucleotides of all ORFs (e.g. from -100 to +50 relative to the start codon of every ORF)', '', 'Usage', '', '  python3 retrieve-seq.py file_name.FASTA file_name.txt [-start/-stop/-startstop rp_from rp_to]', '  [optional commands]', '', 'Input .txt file format example', '', '  ORF-1 13092 13769', '  ORF-2 13738 12953', '  ORF-3 54562 54678', '', "If you don't use the optional commands, the script will retrieve the sequence of every ORF.", 'Using the optional commands you can retrieve a sequence relative to that nt (e.g. -start -100 +150 will retrieve the sequence from -100 to 150 relative to the first nt of every ORF)']
	for item in help_list:
		print(item)
	quit()

def error_msg(stopped = False):
	if stopped: 
		print('retrieve-seq stopped working')
	else: #input error
		print('Wrong input. If you need help, use the -help command (retrieve-seq -help)')
	quit()

def main():
	'''
	Everything starts here.
	'''
	usr_input = sys.argv

	if len(usr_input) < 2:
		error_msg()

	if usr_input[1] in ['-h', '-help']: #retrieve-seq -h/-help
		print_help()

	elif len(usr_input) < 3:
		error_msg()

	#retrieve-seq file_name.FASTA file_name.txt | -start/-stop/-startstop rp_from rp_to
	fasta_fname, plaintext_fname = usr_input[1], usr_input[2]
	try: #for the optional part (relative position)
		rp_option = usr_input[3]
		if rp_option not in ['-start', '-stop', '-startstop']:
			error_msg()
		rp_on = True
	except IndexError:
		rp_on = False

	if rp_on:
		try:
			rp_from, rp_to = int(usr_input[4]), int(usr_input[5])
		except ValueError:
			error_msg()
		except IndexError:
			error_msg()
		retrieve_run(fasta_fname, plaintext_fname, rp_option, rp_from, rp_to)
	elif not rp_on:
		retrieve_run(fasta_fname, plaintext_fname)

	else:
		error_msg()


if __name__ == '__main__':
	main()
