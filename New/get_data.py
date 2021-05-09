import numpy as np
import pandas as pandas
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

filename= "/root/labs/UnderstandingDNA/New/dna_dataset/gbbct1.seq"
def read_genbank(filename):
    '''
    Arg : filepath
    return : Dna sequence
    '''
    sequence = []
    count = 0
    file = open(filename)
    for record in parse(file, "genbank"):
        print(repr(record.seq))
        print("Length of DNA sequence " + str(len(record)))
        coding_dna = record.seq
        sequence.append(coding_dna)
        count += 1
        if count == 5:
            break
            
        
    return coding_dna