#!/usr/bin/env python

########################################################################
#
#	Author: Michele Berselli
#		Harvard Medical School
#		berselli.michele@gmail.com
#
#   Basic objects to work with DNA sequences.
#    The library can use canonical DNA nucleotides and n, N.
#
########################################################################

########################################################################
# Libraries
########################################################################

import sys, os
import math
import zlib

########################################################################
# DNAseq object
#   * most functions are case sensitive *
########################################################################
class DNAseqError(Exception):
    """Custom exception for error tracking for DNAseq.
    """

class DNAseq(object):
    """This is a class to describe a DNA sequence.
    Provide methods to calculate statistics for the sequence.
    """

    def __init__(self, seq):
        """Constructor method.
        Initialize DNAseq object for sequence.

        :param seq: DNA sequence
        :type seq: str

        :ivar seq: DNA sequence [str]
        :ivar len: DNA sequence length [int]
        :ivar ncounts: nucleotide counts [dict]
            case sensitive
        :ivar ncounts_: nucleotide counts [dict]
            case insensitive
        """

        # Attributes
        self.seq = seq
        self.len = len(seq)
        # Raise exception if sequence is null
        if self.len == 0:
            raise DNAseqError('\nSequence is null\n')
        self.ncounts =  None # case sensitive
        self.ncounts_ = None # case insenstive
                             # used for statistics calculation

        # Initialize attributes
        self._count_nucleotides()


    def _count_nucleotides(self):
        """Count number of nucleotides in sequence.
        Include n and N counts.
        """

        # Count per nucleotide, case sensitive
        self.ncounts = {'n': 0, 'n_': 0, 'a': 0, 't': 0, 'c': 0, 'g': 0,
                        'N': 0, 'N_': 0, 'A': 0, 'T': 0, 'C': 0, 'G': 0}

        for n in self.seq:
            self.ncounts[n] += 1

        # Aggregated counts, case insensitive
        self.ncounts_ = {'N': 0, 'N_': 0, 'A': 0, 'T': 0, 'C': 0, 'G': 0}

        # Nucleotides
        self.ncounts_['N'] = self.ncounts['n'] + self.ncounts['N']
        self.ncounts_['A'] = self.ncounts['a'] + self.ncounts['A']
        self.ncounts_['T'] = self.ncounts['t'] + self.ncounts['T']
        self.ncounts_['C'] = self.ncounts['c'] + self.ncounts['C']
        self.ncounts_['G'] = self.ncounts['g'] + self.ncounts['G']


        # Totals excluding n, N
        self.ncounts['n_'] = self.ncounts['a'] + self.ncounts['t'] + self.ncounts['c'] + self.ncounts['g']
        self.ncounts['N_'] = self.ncounts['A'] + self.ncounts['T'] + self.ncounts['C'] + self.ncounts['G']
        self.ncounts_['N_'] = self.ncounts_['A'] + self.ncounts_['T'] + self.ncounts_['C'] + self.ncounts_['G']


    def complement(self):
        """Return sequence complement.
        Case sensitive, n and N are maintained.
        """

        complement = {'n':'n', 'a':'t', 'c':'g', 'g':'c', 't':'a',
                      'N':'N', 'A':'T', 'C':'G', 'G':'C', 'T':'A'}

        return ''.join(complement[n] for n in self.seq)


    def reverse_complement(self):
        """Return sequence reverse complement.
        Case sensitive, n and N are maintained.
        """

        return self.complement()[::-1]


    def hamming_distance(self, seq):
        """Calculate Hamming distance with sequence.
        Case sensitive.

        :param seq: DNA sequence to calculate Hamming distance with
        :type seq: str
        """

        return sum(map(str.__ne__, self.seq, seq))


    def entropy_shannon(self):
        """Calculate Shannon entropy for sequence.
        """

        if not self.ncounts_['N_']: # all sequence is N or n
            return math.nan         # can not calculate results

        # Nucleotides probability
        pA = self.ncounts_['A'] / self.ncounts_['N_']
        pT = self.ncounts_['T'] / self.ncounts_['N_']
        pC = self.ncounts_['C'] / self.ncounts_['N_']
        pG = self.ncounts_['G'] / self.ncounts_['N_']

        # Entropy
        # if 0 not in [pA, pT, pC, pG]:

        # funziona con ogni combinazione di basi
        lA = pA * math.log2(pA)
        lT = pT * math.log2(pT)
        lC = pC * math.log2(pC)
        lG = pG * math.log2(pG)

        return - (lA + lT + lC + lG)

        # return math.nan


    def entropy_topological(self):
        """ """

        return math.nan # placeholder


    def chargaff_fariselli(self):
        """Calculate Chargaff scores (Fariselli) for sequence.
        """

        if not self.ncounts_['N_']: # all sequence is N or n
            return math.nan         # can not calculate results

        # Nucleotides counts
        nA = self.ncounts_['A']
        nT = self.ncounts_['T']
        nC = self.ncounts_['C']
        nG = self.ncounts_['G']

        # Chargaff
        if 0 not in [nA, nT, nC, nG]:
            return abs((nA - nT) / (nA + nT)) + abs((nC - nG) / (nC + nG))

        return math.nan


    def chargaff_taccioli(self):
        """Calculate Chargaff scores (Taccioli) for sequence.
        """

        if not self.ncounts_['N_']: # all sequence is N or n
            return math.nan         # can not calculate results

        # Nucleotides counts
        nA = self.ncounts_['A']
        nT = self.ncounts_['T']
        nC = self.ncounts_['C']
        nG = self.ncounts_['G']

        # Chargaff
        if 0 not in [nA, nT, nC, nG]:
            rAT = min(nA, nT) / max(nA, nT)
            rCG = min(nC, nG) / max(nC, nG)
            return (rAT + rCG) / 2

        return math.nan

    def complexity_kolmogorov(self):
        """Calculate DEFLATE compression.
        That is a variation of LZSS compression algorithm (Lempel–Ziv–Storer–Szymanski).
        """

        seq_ = self.seq.encode()
        return len(zlib.compress(seq_)) / self.len


    def test(self):

        pass
        # r_sh = self.entropy_shannon()
        # r_pf = self.chargaff_fariselli()
        # r_ct = self.chargaff_taccioli()
        # pA = float(self.ncounts_['A'] / self.ncounts_['N_'])
        # pT = float(self.ncounts_['T'] / self.ncounts_['N_'])
        # pC = float(self.ncounts_['C'] / self.ncounts_['N_'])
        # pG = float(self.ncounts_['G'] / self.ncounts_['N_'])
        #
        # return r_pf,r_ct,r_sh,pA,pT,pC,pG


def main():
    """
    """

    pass


if __name__ == '__main__':

    main()

#end if
