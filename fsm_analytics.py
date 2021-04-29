import numpy as np
import pandas as pandas
import get_data
import hmm_viterbi as hmm
import reverse_translator
import finite_state_machine as fsm

filename = '/root/NLU/Project/data/test_cases.fa'
def read_dna(filename):
    '''
    First state of the finite state machine
    '''
    dna_sequence = get_data.read_fasta(filename)
    if dna_sequence == dna_sequence:
        newState = 'Extract Bases' 
    return (newState, dna_sequence)

def extract_base(sequence):
    'Second state of FSM'
    bases = []
    for base in sequence:
        bases.append(base)
    newState = 'Convert Bases'
    return (newState, bases)

def convert_base(bases):
    'Third state of FSM'
    'Transcribes DNA to mRNA'
    cbases = []
    for base in bases:
        if base == 'A':
            cbase = 'U'
        elif base == 'C':
            cbase = 'G'
        elif base == 'G':
            cbase = 'C'
        else:
            cbase = 'A'
        cbases.append(cbase)
    if cbases == cbases:
        newState = 'mRNA'
    return (newState, cbases)

def rna_sequence(cbases):
    'Fourth state in FSM'
    sequence = ''.join([str(cbase) for cbase in cbases])
    newState = 'error_state'
    return (newState, sequence)


m = fsm.FSM()
m.add_state('Read DNA', read_dna, end_state = 0)
m.add_state('Extract Bases', extract_base)
m.add_state('Convert Bases', convert_base)
m.add_state('mRNA', rna_sequence)
m.add_state('error_state', None, end_state = 1)
m.set_start('Read DNA')
m.run(filename)




