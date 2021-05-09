from stateMachine import StateMachine as StateMachine
import numpy as np
import pandas as pandas
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
from get_data import read_genbank as read_genbank


def extract_oneBase(sequence):
    onebase_sequence = []
    coding_dna = list(map(str,sequence))
    #string of DNA sequence 
    sting_dna = str(coding_dna)
    
    #split string to be ['A','ATCG']
    #Append 1 base to list
    onebase_sequence.append(coding_dna[0])
    #Append existing base to list
    ex_sequence = "".join(coding_dna[1:])
    onebase_sequence.append(ex_sequence)
    return onebase_sequence


filename= "/root/labs/UnderstandingDNA/New/dna_dataset/gbbct1.seq"
def read_dna(filename):
    coding_dna = read_genbank(filename)
    if coding_dna == coding_dna:
        newState = "First state"
    return (newState,coding_dna)

def first_state(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'A':
        newState = "state7"
    elif sequence[0] == 'T'or sequence[0]=="C":
        
        newState = "state4"
    elif sequence[0] == "G":
        newState = 'state247'
    return(newState,sequence[1])
    
def state4(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
    
        newState = "state2458"
    elif sequence[0] == 'A'or sequence[0]=="G":
        
        newState = "state58"
    elif sequence[0] == "C":
        
        newState = 'state48'

    return(newState,sequence[1])
    
    
def state7(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state38"
        
    elif sequence[0] == 'A'or sequence[0]=="C":
        
        newState = "state58"
        
    elif sequence[0] == "G":
        
        newState = 'state8'
        
    return(newState,sequence[1])
    
    
def state247(coding_dna):
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state458236"
        
    elif sequence[0] == 'A':
        
        newState = "state36785"
        
    elif sequence[0] == "C":
        
        newState = 'state4872'
    elif sequence[0] == "G":
        newState = 'state5863'
        
    return(newState,sequence[1])
    
def state2458(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        newState = "state458236"
    elif sequence[0] == 'A'or sequence[0]=="G":
        newState = "state5863"
    elif sequence[0] == "C":
        newState = 'state4832'
    return(newState,sequence[1])

def state58(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T' or sequence[0] == 'A'or sequence[0]=="C":
        
        newState = "state38"     
    elif sequence[0] == "G":
        
        newState = 'state8'
    return(newState,sequence[1])

def state48(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state2458"
    elif sequence[0] == 'A':
        
        newState = "state58"
    elif sequence[0]=="C":
        
        newState = "state48"
    elif sequence[0] == "G":
        
        newState = 'state58'
    return(newState,sequence[1])
def state458236(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        newState = "state458236"
    elif sequence[0] == 'A'or sequence[0]=="G":
        newState = "state5863"
    elif sequence[0] == "C":
        newState = 'state4835'
    return(newState,sequence[1])

def state5863(coding_dna):
    sequence = extract_oneBase(coding_dna)
    newState = "Rejection state"
    return(newState,sequence[1])

def state4832(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state45826"
    elif sequence[0] == 'A':
        
        newState = "state5863"
    elif sequence[0]=="C":
        
        newState = "state482"
    elif sequence[0] == "G":
        
        newState = 'state5863'
    return(newState,sequence[1])

def state4872(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state458236"
    elif sequence[0] == 'A':
        
        newState = "state36785"
    elif sequence[0]=="C":
        
        newState = "state4872"
    elif sequence[0] == "G":
        
        newState = 'state5863'
    return(newState,sequence[1])

def state38(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state8"
    elif sequence[0] == 'A':
        
        newState = "state8"
    elif sequence[0]=="C":
        
        newState = "state8"
    elif sequence[0] == "G":
        
        newState = 'state8'
    return(newState,sequence[1])

def state45826(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state458236"
    elif sequence[0] == 'A':
        
        newState = "state5863"
    elif sequence[0]=="C":
        
        newState = "state4832"
    elif sequence[0] == "G":
        
        newState = 'state5863'
    return(newState,sequence[1])

def state482(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state45826"
    elif sequence[0] == 'A':
        
        newState = "state5863"
    elif sequence[0]=="C":
        
        newState = "state482"
    elif sequence[0] == "G":
        
        newState = 'state5863'
    return(newState,sequence[1])

def state378(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state38"
    elif sequence[0] == 'A':
        
        newState = "state78"
    elif sequence[0]=="C":
        
        newState = "state78"
    elif sequence[0] == "G":
        
        newState = 'state8'
    return(newState,sequence[1])

def state36785(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state8"
    elif sequence[0] == 'A':
        
        newState = "state8"
    elif sequence[0]=="C":
        
        newState = "state8"
    elif sequence[0] == "G":
        
        newState = 'state38'
    return(newState,sequence[1])

def state78(coding_dna):
    
    sequence = extract_oneBase(coding_dna)
    print(sequence[0])
    if sequence[0] == 'T':
        
        newState = "state38"
    elif sequence[0] == 'A':
        
        newState = "state78"
    elif sequence[0]=="C":
        
        newState = "state78"
    elif sequence[0] == "G":
        
        newState = 'state8'
    return(newState,sequence[1])

def state8(coding_dna):

    sequence = extract_oneBase(coding_dna)
    newState = "Final state"
    return(newState,sequence[1])



m = StateMachine()

m.add_state("Read DNA",read_dna,end_state=0)
m.add_state("First state",first_state)
m.add_state("state4",state4)
m.add_state("state7",state7)
m.add_state("state8",state8)
m.add_state("state38",state38)
m.add_state("state48",state48)
m.add_state("state58",state58)
m.add_state("state78",state78)
m.add_state("state247",state247)
m.add_state("state378",state378)
m.add_state("state482",state482)
m.add_state("state2458",state2458)
m.add_state("state4872",state4872)
m.add_state("state5863",state5863)
m.add_state("state4832",state4832)
m.add_state("state36785",state36785)
m.add_state("state45826",state45826)
m.add_state("state458236",state458236)
m.add_state('Final state', None, end_state = 1)
m.set_start("Read DNA")
m.run(filename)

'''
m.add_state("Start", start_transitions)
m.add_state("Python_state", python_state_transitions)
m.add_state("is_state", is_state_transitions)
m.add_state("not_state", not_state_transitions)
m.add_state("neg_state", None, end_state=1)
m.add_state("pos_state", None, end_state=1)
m.add_state("error_state", None, end_state=1)
m.set_start("Start")
m.run("Python are great very much a lot")
m.run("Python is difficult")
m.run("Perl is ugly")
'''