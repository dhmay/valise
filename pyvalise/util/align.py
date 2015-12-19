#!/usr/bin/env python
"""Tools for sequence alignment. Quick and dirty, using pairwise2"""

import logging
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from pyvalise.util import dna

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

log = logging.getLogger(__name__)

def global_align_sequences(seq_to_align, ref_seq, matrix=matlist.blosum62,
                           gap_open=-10, gap_extend=-0.5):
    """ Do a pairwise alignment and return the best alignment """
    alns = pairwise2.align.globalds(seq_to_align, ref_seq, matrix,
                                    gap_open, gap_extend)
    return alns[0]


def global_align_sequence_against_sequences(seq_to_align, ref_seqs,
                                            matrix=matlist.blosum62,
                                            gap_open=-10, gap_extend=-0.5):
    return [global_align_sequences(seq_to_align, ref_seq, matrix, gap_open, gap_extend) for ref_seq in ref_seqs]


def find_best_global_alignment(seq_to_align, ref_seqs,
                               matrix=matlist.blosum62,
                               gap_open=-10, gap_extend=-0.5):
    all_alignments = global_align_sequence_against_sequences(seq_to_align, ref_seqs, matrix, gap_open, gap_extend)
    best_alignment = None
    best_score = -1000
    for alignment in all_alignments:
        if alignment[2] > best_score:
            best_alignment = alignment
            best_score = alignment[2]
    return best_alignment


def global_align_sequences_against_sequences(seqs_to_align, ref_seqs,
                                             matrix=matlist.blosum62,
                                             gap_open=-10, gap_extend=-0.5):
    result = list()
    for seq_to_align in seqs_to_align:
        result.append(
            [global_align_sequences(seq_to_align, ref_seq, matrix, gap_open, gap_extend) for ref_seq in ref_seqs])
    return result


def find_different_positions(alignment):
    result = list()
    for i in xrange(0, len(alignment[0])):
        if alignment[0][i] != alignment[1][i]:
            result.append(i)
    return result


def collect_alignments_diffpositions(alignments):
    result = list()
    for alignment in alignments:
        result.extend(find_different_positions(alignment))
    return result


def parse_alignment_string(alignment_string):
    result = list()
    lines = alignment_string.split('\n')
    result.append(lines[0])
    result.append(lines[2])
    if len(lines) > 3:
        score_line = lines[3]
        score = score_line[score_line.find('=') + 1:]
        result.append(score)
    return result


def find_gaps(alignment):
    """find the starting positions of gaps"""
    result = list()
    for i in xrange(0, len(alignment[0])):
        if (alignment[0][i] == '-' or alignment[1][i] == '-') and (i - 1) not in result:
            result.append(i)
    return result


def find_mismatches(alignment):
    result = list()
    for i in xrange(0, len(alignment[0])):
        if alignment[0][i] != '-' and alignment[1][i] != '-' and alignment[0][i] != alignment[1][i]:
            result.append(i)
    return result


def find_mismatch_gap_positions(alignment):
    """find all positions in an alignment that are not exact matches"""
    result = list()
    for i in xrange(0, len(alignment[0])):
        if alignment[0][i] != alignment[1][i]:
            result.append(i)
    return result


def alignment_pretty_string(alignment):
    connector_string_list = list()
    gaps_mismatches = find_mismatch_gap_positions(alignment)
    for i in xrange(0, len(alignment[0])):
        if i in gaps_mismatches:
            connector_string_list.append(' ')
        else:
            connector_string_list.append('|')
    connector_string = ''.join(connector_string_list)
    result = alignment[0] + "\n" + connector_string + "\n" + alignment[1]
    return result


def get_alignment_string(dna_sequence1, dna_sequence2):
    return alignment_pretty_string([dna_sequence1, dna_sequence2])


def write_alignments_to_file(alignments, outfile):
    """write a list of alignments in my own format"""
    i = 0
    for alignment in alignments:
        outfile.write(alignment_pretty_string(alignment) + "\n")
        if i < len(alignments) - 1:
            outfile.write('\n')
        i += 1


def read_alignments_from_file(infile):
    """read a list of alignments in my own format"""
    file_contents = infile.read()
    alignment_chunks = file_contents.split('\n\n')
    alignments = [parse_alignment_string(x) for x in alignment_chunks]
    return alignments


def remove_sequence_gaps(seq):
    return ''.join(seq.split('-'))


def slide_align_perfect(query_seq, ref_seq, match_base_score=1):
    """
    Perform a sliding-only alignment (no gaps) and return the alignment with the longest perfect agreement,
    or None if there is none
    :param query_seq:
    :param ref_seq:
    :return:
    """
    max_perfect_overlap = 0
    max_perfect_overlap_offset = 0
    # ref_seq_offset is the offset between the first position in ref_seq and the first position in query_seq
    for ref_seq_offset in xrange(-len(ref_seq) + 1, len(query_seq) - 1):
        has_mismatch = False
        n_overlap_bases = 0
        for overlap_position in xrange(max(0, ref_seq_offset), min(len(query_seq), ref_seq_offset + len(ref_seq))):
            if query_seq[overlap_position] != ref_seq[overlap_position - ref_seq_offset]:
                has_mismatch = True
                break
            n_overlap_bases += 1
        if not has_mismatch:
            perfect_overlap = n_overlap_bases
            if perfect_overlap > max_perfect_overlap:
                max_perfect_overlap = perfect_overlap
                max_perfect_overlap_offset = ref_seq_offset
    if max_perfect_overlap == 0:
        return None
    ndashes_before_1 = max(0, -max_perfect_overlap_offset)
    ndashes_after_1 = max(0, max_perfect_overlap_offset + len(ref_seq) - len(query_seq))
    alignedseq_1 = ''.join(['-'] * ndashes_before_1) + query_seq + ''.join(['-'] * ndashes_after_1)
    ndashes_before_2 = max(0, max_perfect_overlap_offset)
    ndashes_after_2 = max(0, -max_perfect_overlap_offset + len(query_seq) - len(ref_seq))
    alignedseq_2 = ''.join(['-'] * ndashes_before_2) + ref_seq + ''.join(['-'] * ndashes_after_2)
    return alignedseq_1, alignedseq_2, match_base_score * max_perfect_overlap



    
    
        
        