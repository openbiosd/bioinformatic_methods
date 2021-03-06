import os

#genbank - gbk

def wanted_ids(fasta, txt):
    """
	Given a fasta file and a txt file with wanted sequence ids returns a fasta file with wanted ids.
	"""
    from Bio import SeqIO

    fasta_file = "fasta_file.fasta" # Input fasta file
    wanted_file = "wanted_file.txt" # Input interesting sequence IDs, one per line
    result_file = "result_file.fasta" # Output fasta file

    wanted = set()1
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
	I'm still not sure how it is suppose to work.
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
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline("clustalw2", infile=fasta, pwgapopen=1000, pwgapext=100)
    #help(ClustalwCommandline)
    stdout, stderr = cline()

clustal('dna.fasta')
def allignemt_clustal(aln):
    """
    Takes a fasta file and does a clustal alignment. Requires installed clustal.
    """
    from Bio import AlignIO
    align = AlignIO.read("lab3.aln", "clustal")
    print(align)

def phylo_tree(dnd, draw):
    """
    Takes a dnd file and draws a phyloginetic tree. If draw is False draws ascii tree.
    """
    from Bio import Phylo
    tree = Phylo.read("lab3.dnd", "newick")
    if draw:
        Phylo.draw(tree)
    else:
        Phylo.draw_ascii(tree)

def translate_dna(fasta):
    """
    Translate DNA from fasta to amino acids
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna
    amino_seq = []
    fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
    for seq_record in fasta_sequences:
        seq_rna = seq_record.seq.transcribe()
        seq_dna = seq_rna.translate()
        seq_record.seq = seq_dna
        amino_seq.append(seq_record)
    SeqIO.write(amino_seq, 'amino_seq.fasta', "fasta")
    fasta_sequences.close()

def check_if_dna(fasta):
    """
    Takes a fasta file and checks if sequences are DNA.
    """
    for seq_record in fasta_sequences:
        if all(c.upper() in 'ATGC' for c in seq_record.seq):
            pass # it's DNA

def dna_protein_dna(fasta):
    # Not complete yet. 
    """
    Takes a fasta DNA sequence, converts it to protein, alligns the protein sequence, and converts it back to DNA.
    """
    amino_seq = translate_dna(fasta)
    
