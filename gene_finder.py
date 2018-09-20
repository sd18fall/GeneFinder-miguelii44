# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Miguel Castillo II

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###
def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """  # These 4 tests cover all of the options for nucleotides
    if nucleotide == 'T': #Return opposite codon
        return 'A'
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    genelen = len(dna)
    backward = ""
    while genelen > 0:
        letter = dna[genelen -1]
        complement = get_complement(letter) #Build Complementary Strand
        backward = backward + complement #Add complement to reverse strand
        genelen = genelen-1
    return backward

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    str (dna)
    Strand = ''
    sequence_length = len(dna)
    n = 0
    while n < sequence_length: #Repeat until end of Strand(dna)
        codon = dna[n:n+3] #Breakdown strand into 3's
        if codon == "TAG" or codon == "TAA" or codon == "TGA": #Look for stop codon
            break
        else:
            Strand = Strand + codon
            n = n + 3 #Move forward 3
    return Strand

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    orflist = []
    thrice = 0
    sequence_length = len(dna)
    while thrice < sequence_length:
        start = dna[thrice:thrice+3]
        if start == "ATG": #Look for start codon
            orf = rest_of_ORF(dna[thrice:sequence_length])
            orflist.append(orf) #Add Strand to List
            thrice = thrice + len(orf) #moveforward length of the strand
        else:
            thrice = thrice + 3
    return orflist

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_orf = []
    dnalength = len(dna)
    shiftsingle = dna[1:dnalength]
    shiftdouble = dna[2:dnalength]

    all_orf = find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(shiftsingle) + find_all_ORFs_oneframe(shiftdouble)
    return all_orf

def find_all_ORFs_both_strands(dna):  #Only Function out of order :()
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    #Create a list to hold all ORF strands
    all_orf = []
    reverse = get_reverse_complement(dna)
    revshiftsingle = reverse[1:len(dna)]
    revshiftdouble = reverse[2:len(dna)]
    shiftsingle = dna[1:len(dna)]
    shiftdouble = dna[2:len(dna)]

    all_orf = find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(shiftsingle) + find_all_ORFs_oneframe(shiftdouble) + find_all_ORFs_oneframe(reverse) + find_all_ORFs_oneframe(revshiftsingle) + find_all_ORFs_oneframe(revshiftdouble)
    return all_orf

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_orf = find_all_ORFs_both_strands(dna)
    return max(all_orf)

def load_seq(fasta_file):
    """ Reads a FASTA file and returns the DNA sequence as a string.

    fasta_file: the path to the FASTA file containing the DNA sequence
    returns: the DNA sequence as a string
    """
    retval = ""
    f = open(fasta_file)
    lines = f.readlines()
    for l in lines[1:]:
        retval += l[0:-1]
    f.close()
    return retval

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    for trial in range(num_trials):
        shuffled = shuffle_string(dna)
        all_orf = find_all_ORFs_both_strands(shuffled)
    longest = max(all_orf)
    return len(longest)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    Strand = ''
    n = 0
    while n < len(dna) - 2: #Repeat until end of Strand(dna), cuts off extra letters
        codon = dna[n:n+3] #Breakdown strand into 3's
        amino_acids = aa_table[codon]
        Strand = Strand + amino_acids
        n = n + 3 #Move forward 3
    return Strand


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    both_ORF = find_all_ORFs_both_strands(dna)
    beat = []
    for orf in both_ORF:
        if len(orf) > threshold:
            aminos = coding_strand_to_AA(orf)
            beat.append(aminos)
    return beat
    
dna = load_seq("./data/X73525.fa")
gene_finder(dna)

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
