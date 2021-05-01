import numpy as np
import pandas as pandas
import get_data
import hmm_viterbi as hmm
import reverse_translator
import finite_state_machine as fsm

filename = '/root/NLU/Project/data/test_cases.fa'
def read_dna(filename, pos = 1):
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
    if cbases[0] == 'AUG':
        sequence = cbases
        newState = 'Stop Codon'

    elif cbases[0] == 'A' or cbases[0] == 'U' or cbases[0] == 'C' or cbases[0] == 'G':  
        sequence = ''.join([str(cbase) for cbase in cbases])
        newState = 'Start Codon'
    else:
        newState = 'Final State'
        sequence = cbases
    return (newState, sequence)

def start_codon(sequence):
    ngram = 3
    codons = [sequence[i: i + ngram] for i in range(0, len(sequence), ngram)]
    start = 'AUG'
    start_seq = []
    if codons[0] == start:
        for codon in codons:
            if codon == 'AUG':
                start_seq.append(codon)
        newState = 'Final State'
    else:
        for codon in codons:
            if codon == start:
                start_seq.append(codon)
                codons.remove(codon)
                codons = ['AUG'] + codons
        newState = 'mRNA'
    return (newState, codons)

def stop_codon_state(sequence):
    stop_codons = ['UGA', 'UAA', 'UAG']
    stop_codons_list = []
    for seq in sequence:
        if seq == 'UGA' or seq == 'UAA' or seq == 'UAG':
            stop_codons_list.append(sequence.index(seq))
    newState = 'Final State'
    return (newState, stop_codons_list)


m = fsm.FSM()
m.add_state('Read DNA', read_dna, end_state = 0)
m.add_state('Extract Bases', extract_base)
m.add_state('Convert Bases', convert_base)
m.add_state('mRNA', rna_sequence)
m.add_state('Start Codon', start_codon)
m.add_state('Stop Codon', stop_codon_state)
m.add_state('Final State', None, end_state = 1)
m.set_start('Read DNA')
m.run(filename)



