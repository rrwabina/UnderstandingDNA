from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

file = open("./dna_dataset/gbbct1.seq")
for record in parse(file, "genbank"):
    print(repr(record.seq))
    print("Length of DNA sequence " + str(len(record)))
    coding_dna = record.seq
    #Reverse
    templete_dna = coding_dna.reverse_complement()

    #Transciption to mRNA
    messenger_rna = coding_dna.transcribe()
    print(repr(messenger_rna))

    break



# Create some Grammar Productions
from nltk import CFG
import nltk
grammar = CFG.fromstring("""
S -> Q
Q -> 'A' Q | 'T' Q | 'C' Q | 'G' Q
Q -> P Q
Q -> 
P -> 'A' P 'T' | 'T' P 'A' | 'C' P 'G' | 'G' P 'C'
P -> 
""")

print('A Grammar:', grammar)
print('grammar.start()   =>', grammar.start())
print('grammar.productions() =>')
# Use string.replace(...) is to line-wrap the output.
print(grammar.productions())

sent = ['A', 'T', 'T', 'G', 'A', 'G', 'T']
sent2 = ['C', 'T', 'G', 'C', 'A', 'G', 'C', 'C', 'G', 'C', 'C', 'G', 'A', 'C', 'T', 'G', 'A', 'A', 'A', 'T', 'C', 'T', 'A', 'T', 'C', 'G', 'G', 'G', 'A', 'A', 'G', 'A', 'A', 'A', 'A', 'G', 'C', 'T', 'C', 'G', 'C', 'T', 'T', 'A', 'C', 'G', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'T', 'A', 'A', 'C', 'C', 'C', 'G', 'C', 'A', 'G', 'G', 'A', 'T', 'C', 'C', 'A', 'G', 'T', 'C', 'G', 'C', 'T', 'T', 'A', 'C', 'C', 'T', 'C', 'G', 'C', 'A', 'T', 'C', 'T', 'C', 'A', 'A', 'A', 'A', 'G', 'C', 'A', 'G', 'A', 'A', 'A', 'T', 'A', 'C', 'G', 'G', 'G', 'A', 'G', 'A', 'T', 'A', 'A', 'A', 'C', 'A', 'C', 'A', 'A', 'C', 'T', 'T', 'A', 'T', 'G', 'G', 'T', 'G', 'A', 'G', 'A', 'A', 'C', 'T', 'C', 'C', 'T', 'G', 'T', 'A', 'C', 'C', 'G', 'C', 'T', 'T', 'T', 'A', 'C', 'C', 'T', 'A', 'C', 'G', 'T', 'T', 'G', 'G', 'G', 'C', 'G', 'G', 'T', 'C', 'T', 'C', 'C', 'A', 'T', 'C', 'C', 'T', 'C', 'A', 'G', 'C', 'G', 'T', 'G', 'C', 'T', 'T', 'G', 'C', 'G', 'T', 'T', 'C', 'C', 'T', 'A', 'G', 'C', 'C', 'A', 'T', 'T', 'T', 'G', 'G', 'C', 'A', 'A', 'A', 'T', 'T', 'G', 'C', 'G', 'G', 'C', 'A', 'G', 'C', 'T', 'T', 'C', 'A', 'G', 'G', 'A', 'T', 'T', 'T', 'T', 'T', 'A', 'G', 'G', 'C', 'A', 'A', 'A', 'A', 'C', 'T', 'T', 'T', 'T', 'C', 'C', 'T', 'G', 'G', 'C', 'T', 'C', 'C', 'C', 'T', 'G', 'C', 'G', 'C', 'A', 'C', 'T', 'T', 'T', 'G', 'C', 'A', 'G', 'G', 'A', 'T', 'T', 'T', 'G', 'T', 'T', 'T', 'G', 'G', 'A', 'T', 'G', 'G', 'C', 'T', 'T', 'T', 'C', 'A', 'G', 'A', 'T', 'C', 'C', 'C', 'T', 'T', 'C', 'T', 'T', 'T', 'G', 'A', 'T', 'A', 'A', 'C', 'G', 'G', 'C', 'C', 'C', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T', 'T', 'A',]
parser = nltk.ChartParser(grammar)
for tree in parser.parse(sent2):
     print(tree)