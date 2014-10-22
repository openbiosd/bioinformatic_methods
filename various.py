import os

#help(ClustalwCommandline)

#genbank - gbk

def wanted_ids(fasta, txt):
	"""
	Given a fasta file and a txt file with wanted sequence ids returns a fasta file with wanted ids.
	"""
	from Bio import SeqIO

	fasta_file = "fasta_file.fasta" # Input fasta file
	wanted_file = "wanted_file.txt" # Input interesting sequence IDs, one per line
	result_file = "result_file.fasta" # Output fasta file

	wanted = set()
	with open(wanted_file) as f:
    	for line in f:
    	    line = line.strip()
    	    if line != "":
    	        wanted.add(line)

	fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
	with open(result_file, "w") as f:
    	for seq in fasta_sequences:
    	    if seq.id in wanted:
    	        SeqIO.write([seq], f, "fasta")

def wanted_ids_v2():
	"""                                                                                 
	%prog some.fasta wanted-list.txt                                                    
	"""                                                                                 
	from Bio import SeqIO                                                               
	import sys                                                                          
                                                                                    
	wanted = [line.strip() for line in open(sys.argv[2])]                               
	seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')                                    
	SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.stdout, "fasta")

def reverse_complement(fasta):
	"""
	Takes a fasta file and print the reverse complimantary sequences.
	"""
	from Bio import SeqIO
	for record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    	print record.id
    	print record.seq.reverse_complement()

def print_sequences(fasta):
	"""
	Takes a fasta file and prints all the sequences inside it.
	"""
	from Bio import SeqIO
	for seq_record in SeqIO.parse("lab3.fasta", "fasta"):
		print seq_record.id
		print repr(seq_record.seq)
		print len(seq_record)

def print_ids(fasta):
	"""
	Takes a fasta file and prints all ids in the file.	
	"""
	from Bio import SeqIO
	# used to get aln finles extensions
	handle = open("opuntia.aln", "rU")
	for record in SeqIO.parse(handle, "clustal") :
    	print record.id
	handle.close()

def clustal(fasta):
	"""
	Given a fasta file outputs a MSA as dnd and aln files with same name as fasta.
	"""
	cline = ClustalwCommandline("clustalw2", infile=fasta)
	stdout, stderr = cline()

def allignemt_clustal(aln):
	"""
	Takes a fasta file and does a clustal alignment. Requires installed clustal.
	"""
	from Bio import AlignIO
	align = AlignIO.read("lab3.aln", "clustal")
	print(align)

def phylo_tree(dnd, draw==True):
	"""
	Takes a dnd file and draws a phyloginetic tree. If draw is False draws ascii tree.
	"""
	from Bio import Phylo
	tree = Phylo.read("lab3.dnd", "newick")
	if draw:
		Phylo.draw(tree)
	else:
		Phylo.draw_ascii(tree)
