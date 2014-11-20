#!/usr/bin/env python
""" Utilities for DNA manipulation """

import re
import logging
from Bio.SeqUtils import MeltingTemp


log = logging.getLogger(__name__)

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

# this is copied from http://www.petercollingridge.co.uk/python-bioinformatics-tools/codon-table
DNA_BASES = ['T', 'C', 'A', 'G']
CODONS = [a + b + c for a in DNA_BASES for b in DNA_BASES for c in DNA_BASES]
AMINO_ACIDS_FORTRANSLATE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS_FORTRANSLATE))

# regular expression for checking that a string is DNA. Allow lower and
# upper case
IS_DNA_REX = re.compile("^[ACGTacgt]*$")
IS_DNA_OR_GAP_REX = re.compile("^[ACGTacgt-]*$")
IS_DNA_OR_GAP_OR_N_REX = re.compile("^[ACGNTacgnt-]*$")
IS_DNA_OR_N_REX = re.compile("^[ACGNTacgnt]*$")

def calc_nt_proportion_map(dna_sequences):
    """Build a map from each nucleotide to the proportion of bases in dna_sequences that it represents"""
    # initially holds counts. Converted before return
    result = { "A": 0.0, "T": 0.0, "G": 0.0, "C": 0.0 }

    for seq in dna_sequences:
        for nt in result:
            result[nt] += seq.count(nt)
    grandtotal = float(sum(result.values()))
    for nt in result:
        result[nt] /= grandtotal
    return result

def forward_translate_dna_oneframe(dna_sequence, stop_before_stop=False):
    """forward-translate the DNA sequence in Frame 1.
    If stop_before_stop, end with the last AA before the first in-frame stop codon"""
    verify_is_dna(dna_sequence)
    result = ""
    pos = 0
    while pos <= len(dna_sequence) - 3:
        aa = CODON_TABLE[dna_sequence[pos:pos + 3]]
        if stop_before_stop and aa == '*':
            break
        result = result + CODON_TABLE[dna_sequence[pos:pos + 3]]
        pos += 3
    return result


def forward_translate_dna_threeframes(dna_sequence, stop_before_stop=False):
    """forward-translate the DNA sequence in Frame 1,2,3"""
    result = list([forward_translate_dna_oneframe(dna_sequence, stop_before_stop=stop_before_stop)])
    result.append(forward_translate_dna_oneframe(dna_sequence[1:], stop_before_stop=stop_before_stop))
    result.append(forward_translate_dna_oneframe(dna_sequence[2:], stop_before_stop=stop_before_stop))
    return result


def reverse_complement(dna_sequence):
    return complement(reverse(dna_sequence))


def reverse(dna_sequence):
    """simple string reversal"""
    verify_is_dna(dna_sequence)
    return dna_sequence[::-1]


def complement(dna_sequence):
    verify_is_dna(dna_sequence)
    result_list = list()
    for nt in dna_sequence:
        if nt == '-':
            result_list.append(nt)
        else:
            result_list.append(COMPLEMENT_DICT.get(nt.upper(), ''))

    return "".join(result_list)

def calc_melting_temp_2AT_4GC(dna_sequence):
    """calculate dead-simple melting temperature: 2AT + 4GC.
    No length correction"""
    verify_is_dna(dna_sequence)
    dna_sequence_upper = dna_sequence.upper()
    return 2 * dna_sequence_upper.count('A') + \
           2 * dna_sequence_upper.count('T') + \
           4 * dna_sequence_upper.count('G') + \
           4 * dna_sequence_upper.count('C')


def calc_melting_temp_quikchange(dna_sequence):
    """use the Quikchange algorithm to calc melting temperature of a
    fully-matched sequence"""
    return calc_melting_temp_with_mismatches(dna_sequence, dna_sequence)


def calc_melting_temp_with_mismatches(dna_sequence1, dna_sequence2):
    """From the Stratagene Quikchange kit manual here:
    http://www.chem.uky.edu/courses/che554/quikchange.pdf
    formula from the Stratagene
    Quikchange kit: Tm = 81.5 + 0.41(%GC) - 675/N - %mismatch
    
    * N is the primer length in bases
    * %GC is based on matching bases only
    * values for %GC and %mismatch are whole numbers -- so 75% is 75, not .75
    
    Yes, this is super weird. Stratagene notes not appropriate for
    primers longer than 45 bases; plus, it's really for mutagenesis 
    primer design. Need to revisit.
    
    Do not use this on very short sequences, either! 3 is definitely too short.
    """
    N = len(dna_sequence1)
    assert (N == len(dna_sequence2))
    n_matches = 0
    n_matched_gcs = 0

    for i in xrange(0, N):
        if dna_sequence1[i] == dna_sequence2[i]:
            n_matches += 1
            if dna_sequence1[i] == 'G' or dna_sequence2[i] == 'C':
                n_matched_gcs += 1
    # corner case
    if n_matches == 0:
        return 0
    percent_gc = int(100 * float(n_matched_gcs) / float(n_matches))
    percent_mismatch = int(100 * float(N - n_matches) / float(N))
    tm = 81.5 + 0.41 * percent_gc - float(675) / float(N) - percent_mismatch
    log.debug('\n' + get_alignment_string(dna_sequence1, dna_sequence2))
    log.debug(
        "calc_melting_temp_with_mismatches: N=%d, n_matches=%d, n_matched_gcs=%d, percent_gc=%d, percent_mismatch=%d. Tm=%f" %
        (N, n_matches, n_matched_gcs, percent_gc, percent_mismatch, tm))
    return tm


def get_alignment_string(dna_sequence1, dna_sequence2):
    """Create an alignment string between two sequences of equal length, assumed lined up optimally"""
    assert (len(dna_sequence1) == len(dna_sequence2))
    connector_string_list = list()
    for i in xrange(0, len(dna_sequence1)):
        if dna_sequence1[i] == dna_sequence2[i]:
            connector_string_list.append('|')
        else:
            connector_string_list.append(' ')
    connector_string = ''.join(connector_string_list)
    return dna_sequence1 + "\n" + connector_string + "\n" + dna_sequence2


def calc_melting_temp(dna_sequence, dnac=50, saltc=50):
    """from Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594.
    dnac is DNA concentration [nM]
    saltc is salt concentration [mM]."""
    return MeltingTemp.Tm_staluc(dna_sequence, dnac, saltc)


def show_dna_protein_aligned(dna_sequence, protein_sequence):
    """does a simple check to make sure dna_sequence codes for protein_sequence
    and then shows the two aligned. Return the two-line alignment string"""
    verify_is_dna(dna_sequence)
    if not forward_translate_dna_oneframe(dna_sequence) == protein_sequence:
        raise ValueError('DNA sequence does not code for protein sequence:\n' +
                         dna_sequence + "\n" + protein_sequence +
                         "\nInstead, codes for:\n" +
                         forward_translate_dna_oneframe(dna_sequence))
    return '  '.join(protein_sequence) + "\n" + dna_sequence


def reverse_trans_codon_choices_from_table(revtrans_table, protein_seq, codon_choices):
    """reverse-translate a set of codon choices from a reverse translation table"""
    #    print(str(i) for i in range(0,len(protein_seq)))
    return ''.join(revtrans_table[protein_seq[i]][codon_choices[i]] for i in range(0, len(protein_seq)))


def verify_is_dna(input_dna, allow_dashes=True, allow_N=True):
    """Check that an input string is DNA. If so, do nothing. not, raise a
    friendly ValueError"""
    rex = IS_DNA_REX
    if allow_dashes and allow_N:
        rex = IS_DNA_OR_GAP_OR_N_REX
    else:
        if allow_dashes:
            rex = IS_DNA_OR_GAP_REX
        else:
            if allow_N:
                rex = IS_DNA_OR_N_REX
    m = rex.match(input_dna)
    if not m:
        raise ValueError("Input sequence is not DNA. Sequence: %s" % input_dna)


def hamdist(str1, str2, max_care_about=1000000):
    """hamming distance calc, with a maximum that we care about. If we get to
    that maximum, return it (to save compute time)"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
            if diffs >= max_care_about:
                break
    return diffs


