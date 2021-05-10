"RNA_fsm"
#%%
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import numpy as np
import pandas as pandas
import get_data
import hmm_viterbi as hmm
import finite_state_machine as fsm
print('import OK')

#%%
filename = '../data/test_cases.fa'
def read_dna(filename, pos = 0):
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
    print(bases)
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
        newState = 's1'
    return (newState, cbases)
    
def state1(cbases):
    print(cbases[0])
    if cbases[0] == 'U':
        newState = 's2'
    else:
        newState = 's1'
    return (newState, cbases[1:])

def state2(cbases):
    print(cbases[0])
    if cbases[0] == 'A':
        newState = 's3'
    elif cbases[0] == 'G':
        newState = 's4'
    else:
        newState = 's1'
    return (newState, cbases[1:])

def state3(cbases):
    print(cbases[0])
    if cbases[0] == 'A' or cbases[0] == 'G':
        newState = 's5'
    else:
        newState = 's1'
    return (newState, cbases[1:])

def state4(cbases):
    print(cbases[0])
    if cbases[0] == 'A':
        newState = 's5'
    else:
        newState = 's1'
    return (newState, cbases[1:])

def state5(cbases):
    newState = "Final state"
    return (newState, cbases[1:])

#%%
m = fsm.FSM()
m.add_state('Read DNA', read_dna, end_state = 0)
m.add_state('Extract Bases', extract_base)
m.add_state('Convert Bases', convert_base)
m.add_state('s1', state1)
m.add_state('s2', state2)
m.add_state('s3', state3)
m.add_state('s4', state4)
m.add_state('s5', state5)
m.add_state('Final state', None, end_state = 1)
m.set_start('Read DNA')

m.run(filename)
# %%
