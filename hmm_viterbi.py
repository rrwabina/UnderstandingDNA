import numpy as np
import random
import math
import get_data

def prior_nucleotide(sequence, molecule = 'dna'):
    total_seq = len(sequence)
    total_nucleotides = {}
    nucleotides = ['A', 'T', 'C', 'G']
    if molecule == 'dna':
        bases = ['A', 'T', 'C', 'G']
        for base, nuclei in zip(bases, nucleotides):
            total_nucleotides[nuclei] = sequence.count(base)
    elif molecule == 'rna':
        bases = ['A', 'T', 'C', 'G']
        for base, nuclei in zip(bases, nucleotides):
            total_nucleotides[nuclei] = sequence.count(base)
    else:
        return 'Wrong input!'
    low = {"A": 0.291, "T": 0.291, "C": 0.209, "G": 0.209}
    high = {}
    emission = []
    for nuclei in nucleotides:
        high[nuclei] = np.round(total_nucleotides[nuclei]/total_seq, 4)
    emission.append(low)
    emission.append(high)
    return emission

class HMM():
    
    def __init__(self, sequence):
        self.sequence = sequence
        self.num_states = 2
        self.transition = [[0.5, 0.5], [0.5, 0.5]]
        self.prior = [0.5, 0.5]
        self.emission = prior_nucleotide(sequence, molecule = 'dna')
    
    def logprob(self, sequence, states):
        result = []
        initial = self.emission[states[0]][sequence[0]]
        result.append(math.log(self.prior[0]) + math.log(initial))

        for i in range(1, len(sequence)):
            em = self.emission[states[i]][sequence[i]]
            trns = self.transition[states[i-1]][states[i]]
            prev = result[i-1]

            result.append(math.log(trns) + math.log(em) + prev)

        return result[len(result)-1]

    def viterbi(self, sequence):
        l = len(sequence)

        initial0 = self.emission[0][sequence[0]]
        initial1 = self.emission[1][sequence[0]]

        row0 = np.empty([self.num_states, l])
        row1 = np.empty([self.num_states, l])

        row0[0, 0] = math.log(self.prior[0]) + math.log(initial0)
        row0[1, 0] = math.log(self.prior[1]) + math.log(initial1)
        row1[0, 0] = 0
        row1[1, 0] = 0

        for i in range(1, l):
            for j in range(0, self.num_states):

                prev1 = row0[0][i-1]
                prev2 = row0[1][i-1]

                emission = self.emission[j][sequence[i]]

                trns1 = self.transition[0][j]
                trns2 = self.transition[1][j]

                low = math.log(trns1) + prev1
                high = math.log(trns2) + prev2
                row0[j, i] = max(low, high) + math.log(emission)
                row1[j, i] = np.argmax([low, high])

        states = np.empty(l, int)
        states[l-1] = row0[:, l-1].argmax()

        for j in range(l-1, 0, -1):
            states[j-1] = row1[states[j], j]

        return states.tolist()


