#!/usr/bin/env python
"""
Analysis of HTS data is done based on a .bam file. Aligned and unaligned reads
are analyzed separately, but with the same goal: create a dataset with an
entry about each unique coding sequence observed, with information about
what we observed.
"""

import logging
from pyvalise.io import fasta_io
from os.path import basename
from pyvalise.util import dna
import pysam
from pyvalise.util import charts
import csv
import math


__author__ = "Damon May"
__copyright__ = "Copyright (c) Damon May" 
__license__ = ""
__version__ = ""

# Defaults

# minimum worst base quality to keep
DEFAULT_MIN_WORST_BASE_QUAL = 50

# maximum number of differences in an alignment to keep
DEFAULT_MAX_ALIGN_DIFFS = 4
# read length
DEFAULT_READ_LENGTH = 149

# minimum number of reads observed for us to keep a unique sequence
DEFAULT_MIN_READS_TO_KEEP_IN_ANALYSIS = 2

# minimum length of a coded protein to keep in output
MIN_OUTPUT_PROTEIN_LENGTH = 7

DEFAULT_N_CORES = 6

log = logging.getLogger(__name__)

# qnames of all reads
all_read_qnames = set()


def filter_samfile_min_coding_base_qual(in_samfile, out_filepath,
                                        fasta_file,
                                        fiveprime_codingstartseq='GGATCC',
                                        threeprime_aftercodingseq='TAATGC',
                                        min_qualscore=DEFAULT_MIN_WORST_BASE_QUAL):
    """Filter a .bam file based on the *lowest* base quality within the
    coding region, starting with fiveprime_codingstartseq and ending with
    the first stop codon. If no coding sequence, toss it. Also toss too-short
    coding sequences. NOTE: this approach assumes that we know the length of the coding sequence
    (i.e., we're dealing with a variant library). If not, we'll need to determine some other way, e.g.,
    stop codon in the reference sequence."""
    chrom_name_seq_map = fasta_io.load_fasta_name_sequence_map(fasta_file)

    out_samfile = pysam.Samfile(out_filepath, 'wb', template=in_samfile)
    n_reads_kept = 0
    n_reads_evaluated = 0
    n_rej_unmapped = 0
    n_rej_qualscore = 0
    n_rej_nocodingseq = 0
    for aread in in_samfile:
        n_reads_evaluated += 1
        if n_reads_evaluated % 1000000 == 0:
            print("Evaluated %d reads for filter..." % n_reads_evaluated)
        if aread.is_unmapped:
            n_rej_unmapped += 1
            continue
        try:
            chromname = in_samfile.getrname(aread.tid)
            start_codingseq_ind = aread.seq.index(fiveprime_codingstartseq)
            start_codingseq_refpos = aread.positions[start_codingseq_ind]
            end_codingseq_refpos = chrom_name_seq_map[chromname].index(threeprime_aftercodingseq)
            aftercodingseq_ind = 0
            for i in xrange(start_codingseq_ind, len(aread.positions)):
                if aread.positions[i] <= end_codingseq_refpos:
                    aftercodingseq_ind = i
                else:
                    break
        except ValueError:
            n_rej_nocodingseq += 1
            continue
        if aftercodingseq_ind - start_codingseq_ind < 15:
            n_rej_nocodingseq += 1
            continue

        codingseq_qualscores = [ord(x) for x in aread.qual][start_codingseq_ind:aftercodingseq_ind]
        if min(codingseq_qualscores) >= min_qualscore:
            n_reads_kept += 1
            out_samfile.write(aread)
        else:
            n_rej_qualscore += 1
    print("kept %d of %d reads" % (n_reads_kept, n_reads_evaluated))
    print("Rejection reasons: No coding seq = %d, unmapped = %d, low score = %d" %
          (n_rej_nocodingseq, n_rej_unmapped, n_rej_qualscore))
    in_samfile.close()
    out_samfile.close()


def analyze_worst_coding_quality_scores(samfile, fasta_file,
                                        fiveprime_codingstartseq='GGATCC',
                                        threeprime_aftercodingseq='TAATGC'):
    """Return the distributions of *worst* positional quality scores for sequences
    that match the target fasta file, vs. sequences that don't. """
    fasta_seqs = fasta_io.load_fasta_name_sequence_map(fasta_file).values()
    fasta_codingseqs = set(
        [calc_codingseq_safe(seq, fiveprime_codingstartseq, threeprime_aftercodingseq) for seq in fasta_seqs])

    perfect_read_lowestscores = list()
    imperfect_read_lowestscores = list()
    for aread in samfile:
        pos_qualscores = [ord(x) for x in aread.qual]
        read_codingseq = calc_codingseq_safe(aread.seq,
                                             fiveprime_codingstartseq,
                                             threeprime_aftercodingseq)
        if not read_codingseq:
            imperfect_read_lowestscores.append(min(pos_qualscores))
            continue
        lowest_codingqual = min(
            pos_qualscores[aread.seq.index(fiveprime_codingstartseq):aread.seq.index(threeprime_aftercodingseq)])
        if read_codingseq in fasta_codingseqs:
            perfect_read_lowestscores.append(lowest_codingqual)
        else:
            imperfect_read_lowestscores.append(lowest_codingqual)
    return (perfect_read_lowestscores, imperfect_read_lowestscores)


def analyze_htseq_library_run(bam_filename,
                              fiveprime_beforecodingseq,
                              threeprime_aftercodingseq,
                              fasta_file,
                              min_reads_to_keep=DEFAULT_MIN_READS_TO_KEEP_IN_ANALYSIS,
                              remove_outofframe=False,
                              min_worst_base_qual=DEFAULT_MIN_WORST_BASE_QUAL,
                              maxreads=999999999999,
                              should_include_unmapped=True):
    """
    Analyze a HTS run of a scaffold-based library sample, returning an HTSLibraryAnalysis of the unique observed
    sequences that align to members of the provided library fasta (fasta_file). Optionally include unmapped reads.

    NOTE: relies on 5' coding DNA and 3' DNA immediately after coding sequence to be the same for all members.

    :param bam_filename:
    :param fiveprime_codingstartseq:
    :param threeprime_aftercodingseq:
    :param fasta_file: sequences are DNA sequences that are supposed to be in the sample
    :param min_reads_to_keep: minimum reads for a particular sequence in order to keep it
    :param remove_outofframe: remove out-of-frame sequences? (boolean)
    :param min_worst_base_qual:
    :param maxreads: Maximum number of reads to analyze. This is for debugging -- default is a huge number.
    :param should_include_unmapped: (boolean) Handle both mapped and unmapped reads. Unmapped reads are handled in a
    separate method, but the data are combined.
    :return: an HTSLibraryAnalysis object
    """

    # map from "chromosome" names (fasta names) to reference sequences
    chrom_name_seq_map = fasta_io.load_fasta_name_sequence_map(fasta_file)

    # index of the coding sequence end depends on coding sequence length, must be ascertained for every chrom
    # dhmay adding this to address chip 7 variable-length issue
    chromname_codingseqend_map = {}
    for chromname in chrom_name_seq_map:
        chromname_codingseqend_map[chromname] = len(chrom_name_seq_map[chromname])

    samfile = pysam.Samfile(bam_filename, "rb")

    #result raw materials
    seq_obscodingseq_map = dict()

    #accounting
    n_neg_reads = 0
    n_pos_reads = 0
    n_kept_neg_reads = 0
    n_kept_pos_reads = 0
    n_fullcover_reads = 0
    n_reads_checked = 0
    n_reads_mapped = 0
    n_reads_unmapped = 0
    n_unmapped_sequences = 0
    n_reads_kept_mapped = 0
    n_reads_kept_unmapped = 0
    n_rechar_unmapped_mapped = 0

    # this is a strange bit of syntax: if until_eof is true, unmapped reads will be included. If not, it will
    # stop after mapped reads. That took me a while to figure out!
    for aread in samfile.fetch(until_eof=should_include_unmapped):
        if n_reads_checked >= maxreads:
            print("WARNING! Stopping because analyzed %d reads!" % n_reads_checked)
            break
        if n_reads_checked % 50000 == 0:
            print(
            "Checked %d Reads. Mapped: %d. Unmapped: %d. Kept mapped: %d. Kept unmapped: %d. Full-cover: %d. seqs: %d. Unmapped seqs: %d. Rechar unmap as map: %d" %
            (n_reads_checked, n_reads_mapped, n_reads_unmapped,
             n_reads_kept_mapped, n_reads_kept_unmapped,
             n_fullcover_reads, len(seq_obscodingseq_map),
             n_unmapped_sequences, n_rechar_unmapped_mapped))
        #this isn't quite right because of the ambiguity of unmapped reads,
        #but it's just for status message, so who cares?
        if aread.is_reverse:
            n_neg_reads += 1
        else:
            n_pos_reads += 1
        read_qualscores = [ord(x) for x in aread.qual]

        n_reads_checked += 1
        #Handle mapped and unmapped separately
        if aread.is_unmapped:
            # unmapped
            n_reads_unmapped += 1
            if should_include_unmapped:
                read_codingseq = analyze_unmapped_read(aread,
                                                       fiveprime_beforecodingseq,
                                                       threeprime_aftercodingseq,
                                                       min_worst_base_qual=min_worst_base_qual)
                if not read_codingseq:
                    continue
                n_fullcover_reads += 1
            else:
                continue
        else:
            # Mapped read.
            # do a bunch of checks to see if we keep this

            # dhmay changing to address variable lengths
            chromname = samfile.getrname(aread.tid)
            chrom_codingseq_end = chromname_codingseqend_map[chromname]

            n_reads_mapped += 1
            if aread.pos > 0 or aread.aend < chrom_codingseq_end:
                #print("%d  %d, %d" % (aread.pos, aread.aend, chrom_codingseq_end))
                # read starts after coding sequence or ends before
                continue
            # determine the start and end positions of the read relative to the coding sequence. If either of those
            # positions is missing, throw this read away.
            read_coding_start = -1
            read_coding_end = -1
            for i in xrange(0, len(aread.positions)):
                if aread.positions[i] == 0:
                    read_coding_start = i
                if aread.positions[i] == chrom_codingseq_end-1:
                    read_coding_end = i

            if read_coding_start == -1 or read_coding_end == -1:
                # start or end position is missing
                continue
            # check quality scores of all positions in the coding sequence
            lowest_codingqual = min(read_qualscores[read_coding_start:read_coding_end + 1])
            if lowest_codingqual < min_worst_base_qual:
                continue
            # if we got here, read has positions mapping to start and end
            # of coding sequence. There may be gaps
            n_fullcover_reads += 1

            read_codingseq = aread.seq[read_coding_start:read_coding_end + 1]

        if not read_codingseq:
            # this is an unmapped read that we couldn't characterize
            continue

        # get rid of ambiguous sequences.
        if 'N' in read_codingseq or '_' in read_codingseq:
            continue

#        if len(chrom_codingseq) % 3 != len(read_codingseq) % 3:
#            #                    print("OUT OF FRAME")
#            if remove_outofframe:
#                continue

        if not aread.is_unmapped and not read_codingseq:
            print("AARGS: %s %d %d" % (aread.seq, read_coding_start, read_coding_end))

        # OK, this is one to keep

        # update accounting variables
        if aread.is_reverse:
            n_kept_neg_reads += 1
        else:
            n_kept_pos_reads += 1
        if aread.is_unmapped:
            n_reads_kept_unmapped += 1
        else:
            n_reads_kept_mapped += 1
        # accounting. If we're going to reannotate an unmapped sequence as mapped, note the change
        if not aread.is_unmapped and read_codingseq in seq_obscodingseq_map and \
                not seq_obscodingseq_map[read_codingseq].has_refseq:
            n_rechar_unmapped_mapped += 1
            n_unmapped_sequences -= 1
        if not read_codingseq in seq_obscodingseq_map:
            # put a dummy entry in the dict. This entry is appropriate
            # if the read is unmapped
            obscodingseq = ObservedCodingSeq(read_codingseq, 0, 0,
                                             False, '', '', [])
            seq_obscodingseq_map[read_codingseq] = obscodingseq
            if aread.is_unmapped:
                n_unmapped_sequences += 1
        obscodingseq = seq_obscodingseq_map[read_codingseq]

        seq_obscodingseq_map[read_codingseq].readcount += 1

        # now, if the read is mapped and the obscodingseq is not (i.e.,
        # either it was built based on a previous unmapped read, or 
        # it's new for this read), add information about the alignment
        if not aread.is_unmapped and not obscodingseq.has_refseq:
            obscodingseq.has_refseq = True
            obscodingseq.refseq_name = samfile.getrname(aread.tid)
            obscodingseq.refseq_sequence = chrom_name_seq_map[obscodingseq.refseq_name]
            ref_sequence = obscodingseq.refseq_sequence

            obscodingseq.mut_positions_refbased = list()
            obscodingseq.ins_positions_after_refbased = list()
            obscodingseq.del_positions_refbased = list()

            refpos_readpos_map = {}
            last_ref_pos = -1
            for (read_pos, ref_pos) in aread.aligned_pairs:
                if not ref_pos:
                    if last_ref_pos not in obscodingseq.ins_positions_after_refbased:
                        obscodingseq.ins_positions_after_refbased .append(read_pos)
                else:
                    last_ref_pos = ref_pos
                refpos_readpos_map[ref_pos] = read_pos
                if read_pos is not None and ref_pos is not None and aread.seq[read_pos] != ref_sequence[ref_pos]:
                    obscodingseq.mut_positions_refbased .append(ref_pos)

            for refpos in xrange(0, len(ref_sequence)):
                if not refpos in refpos_readpos_map:
                    obscodingseq.del_positions_refbased.append(refpos)

            all_bad_poses_set = set()

            all_bad_poses_set.update(obscodingseq.mut_positions_refbased)
            all_bad_poses_set.update(obscodingseq.ins_positions_after_refbased)
            all_bad_poses_set.update(obscodingseq.del_positions_refbased)
            obscodingseq.all_bad_positions_refbased = list(all_bad_poses_set)
            obscodingseq.all_bad_positions_refbased.sort()

        if not aread.is_reverse:
            obscodingseq.readcount_posstrand += 1

    log.debug("Reads: %d, Full seq: %d" %
              (n_reads_checked, n_fullcover_reads))
    log.debug("Strand: kept %d of %d pos. Kept %d of %d neg" %
              (n_kept_pos_reads, n_pos_reads, n_kept_neg_reads, n_neg_reads))

    return HTSLibraryAnalysis(seq_obscodingseq_map.values(),
                              min_reads=min_reads_to_keep)


def analyze_unmapped_read(aread,
                          fiveprime_beforecodingseq,
                          threeprime_aftercodingseq,
                          min_worst_base_qual=DEFAULT_MIN_WORST_BASE_QUAL):
    """ process an unmapped read by inferring the coding sequence. Very different procedure than for mapped
    reads. If we can do it, return the coding sequence and update the is_reverse
    flag on the read. Otherwise return None"""
    if '_' in aread.seq or 'N' in aread.seq:
        return None
    # note: this assumes fiveprime_beforecodingseq is a palindrome
    if not fiveprime_beforecodingseq in aread.seq:
        return None
    read_revcomp = dna.reverse_complement(aread.seq)
    possible_seqs_with_codingseq = list()
    for read_seq in [aread.seq, read_revcomp]:
        # march back through the sequence that should precede the coding
        # start and look for disagreements
        coding_startpos = read_seq.index(fiveprime_beforecodingseq)
        found_conflict = False
        for i in xrange(0, len(fiveprime_beforecodingseq)):
            pos_to_consider = coding_startpos - i - 1
            if pos_to_consider > 0 and read_seq[pos_to_consider] != fiveprime_beforecodingseq[
                                len(fiveprime_beforecodingseq) - i - 1]:
                found_conflict = True
                break
        if not found_conflict:
            possible_seqs_with_codingseq.append(read_seq)
    if len(possible_seqs_with_codingseq) != 1:
        return None
    seq_with_codingseq = possible_seqs_with_codingseq[0]
    coding_startpos = seq_with_codingseq.index(fiveprime_beforecodingseq) + len(fiveprime_beforecodingseq)
    read_is_reverse_strand = seq_with_codingseq == read_revcomp
    stop_codon_index = coding_startpos
    found_stop = False

    while len(seq_with_codingseq) > stop_codon_index + 3:
        codon = seq_with_codingseq[stop_codon_index:stop_codon_index + 3]
        if codon == threeprime_aftercodingseq[0:3]:
            found_stop = True
            break
        stop_codon_index += 3
    if not found_stop:
        return None
    read_qualscores = [ord(x) for x in aread.qual]
    if min(read_qualscores[coding_startpos:stop_codon_index]) < min_worst_base_qual:
        return None

    #OK, it's a keeper
    read_codingseq = seq_with_codingseq[coding_startpos:stop_codon_index]
    aread.is_reverse = read_is_reverse_strand

    return read_codingseq


def calc_codingseq_safe(dnaseq,
                        fiveprime_codingstartseq='GGATCC',
                        threeprime_aftercodingseq='TAATGCGGCCGC'):
    """
    Calculate the coding sequence given a DNA sequence

    :param dnaseq:
    :param fiveprime_codingstartseq:
    :param threeprime_aftercodingseq:
    :return: Coding sequence, or None if fiveprime coding start sequence *and* 3' sequence after codingsequence
    aren't both found
    """
    if not fiveprime_codingstartseq in dnaseq or \
            not threeprime_aftercodingseq in dnaseq:
        return None
    codingseq_start = dnaseq.index(fiveprime_codingstartseq)
    codingseq_end = dnaseq.index(threeprime_aftercodingseq) - 1
    if codingseq_end <= codingseq_start:
        return None
    return dnaseq[codingseq_start:codingseq_end + 1]


def calc_codingseq(dnaseq,
                   fiveprime_codingstartseq='GGATCC',
                   threeprime_aftercodingseq='TAATGCGGCCGC'):
    """
    Calculate the coding part of a DNA sequence. Assertions will fail
    if doesn't include ends

    :param dnaseq:
    :param fiveprime_codingstartseq:
    :param threeprime_aftercodingseq:
    :return:
    """
    assert (fiveprime_codingstartseq in dnaseq and
            threeprime_aftercodingseq in dnaseq)
    codingseq_start = dnaseq.index(fiveprime_codingstartseq)
    codingseq_end = dnaseq.index(threeprime_aftercodingseq) - 1
    assert (codingseq_end > codingseq_start)
    return dnaseq[codingseq_start:codingseq_end + 1]


def build_unmapped_codingseqs_countmap(samfile, n_reads_to_sample=500000,
                                       fiveprime_codingstartseq='GGATCC',
                                       threeprime_aftercodingseq='TAATGCGGCCGC'):
    """build a map from *unmapped* coding sequences (bounded by known ends)
    to counts of reads. This is for casual use"""
    result = dict()
    n_no_codingseq = 0
    for _ in xrange(0, n_reads_to_sample):
        aread = samfile.next()
        seq = aread.seq
        if not aread.is_unmapped:
            continue
        if fiveprime_codingstartseq in seq and threeprime_aftercodingseq in seq and \
                                seq.index(fiveprime_codingstartseq) + 30 < seq.index(threeprime_aftercodingseq):
            codingseq = seq[seq.index(fiveprime_codingstartseq):seq.index(threeprime_aftercodingseq)]
            if not codingseq in result:
                result[codingseq] = 0
            result[codingseq] += 1
        else:
            n_no_codingseq += 1
    print("No coding sequence: %d of %d" % (n_no_codingseq, n_reads_to_sample))
    return result


def infer_fastq_proteins(fastq_filepath,
                         coding_startseq='GGATCC',
                         n_reads_to_sample=999999999999999,
                         min_reads=1,
                         min_protein_len=1):
    """
    Collect unique protein sequences coded for by fastq DNA, by assuming the coding sequence starts with GGATCC.
    This is hokey, just grabs n_reads_to_sample from the top (huge number by default, so the whole file). Assumes
    + strand.

    Note: inelegant code duplication with infer_fastq_coding_seqs

    :param fastq_filepath:
    :param coding_startseq:
    :param n_reads_to_sample: Number of reads from the top to sample
    :param min_reads:
    :param min_protein_len:
    :return:
    """
    fastq = pysam.Fastqfile(fastq_filepath)
    seq_readcount_map = dict()
    n_had_startseq = 0

    n_sampled = 0
    while True:
        if n_sampled >= n_reads_to_sample:
            break
        try:
            aread = fastq.next()
        except Exception:
            break
        n_sampled += 1
        if not coding_startseq in aread.sequence:
            continue
        n_had_startseq += 1
        seq_coding_to_end = aread.sequence[aread.sequence.index(coding_startseq):]
        try:
            protein_seq = dna.forward_translate_dna_firstframe(seq_coding_to_end,
                                                             stop_before_stop=True)
        except Exception:
            continue
        if len(protein_seq) < min_protein_len:
            continue
        if not protein_seq in seq_readcount_map:
            seq_readcount_map[protein_seq] = 0
        seq_readcount_map[protein_seq] += 1

    print("%d of %d reads had %s" % (n_had_startseq, n_sampled, coding_startseq ))
    result = set()
    for protein_seq in seq_readcount_map:
        if seq_readcount_map[protein_seq] >= min_reads:
            result.add(protein_seq)

    return result


def infer_fastq_coding_seqs(fastq_filepath,
                            coding_startseq='GGATCC',
                            noncoding_startseq='TAATGCGGCCGC',
                            min_codingseq_length=50,
                            n_reads_to_sample=999999999999999,
                            min_reads=1,
                            full_only=True):
    """Infer the coding sequences represented by fastq reads.
    This is hokey, just grabs n_reads_to_sample from the top (default to a huge number, i.e., whole file).

     Note: inelegant code duplication with infer_fastq_proteins"""
    fastq = pysam.Fastqfile(fastq_filepath)
    seq_readcount_map = dict()
    n_full = 0

    n_sampled = 0
    while True:
        if n_sampled >= n_reads_to_sample:
            break
        try:
            aread = fastq.next()
        except Exception:
            break
        n_sampled += 1
        if coding_startseq in aread.sequence:
            if noncoding_startseq in aread.sequence:
                if aread.sequence.index(noncoding_startseq) - aread.sequence.index(
                        coding_startseq) >= min_codingseq_length:
                    n_full += 1
                    codingseq = aread.sequence[
                                aread.sequence.index(coding_startseq):aread.sequence.index(noncoding_startseq)]
                    # remove reads with ambiguous positions
                    if not 'N' in codingseq:
                        if not codingseq in seq_readcount_map:
                            seq_readcount_map[codingseq] = 0
                        seq_readcount_map[codingseq] += 1
            elif not full_only:
                codingseq = aread.sequence[aread.sequence.index(coding_startseq):]
                if len(codingseq) > min_codingseq_length and not 'N' in codingseq:
                    if not codingseq in seq_readcount_map:
                        seq_readcount_map[codingseq] = 0
                    seq_readcount_map[codingseq] += 1
        elif not full_only and noncoding_startseq in aread.sequence:
            codingseq = aread.sequence[:aread.sequence.index(noncoding_startseq)]
            if len(codingseq) > min_codingseq_length and not 'N' in codingseq:
                if not codingseq in seq_readcount_map:
                    seq_readcount_map[codingseq] = 0
                seq_readcount_map[codingseq] += 1
    print("%d of %d had full sequence" % (n_full, n_sampled))
    result = set()
    for seq in seq_readcount_map:
        if seq_readcount_map[seq] >= min_reads:
            result.add(seq)
    return result


def infer_fastq_noncoding_starts_ends(fastq_filepath,
                                      coding_startseq='GGATCC',
                                      noncoding_startseq='TAATGCGGCCGC',
                                      min_end_len=10,
                                      max_end_len=70,
                                      n_reads_to_sample=999999999999,
                                      codingseqs_for_inference=None,
                                      min_longest_support=100,
                                      strand='either'):
    """
    Pull out the most common 5' and 3' sequences outside the coding
    sequences in a fastq file full of reads. The aim is to infer what
    vector the coding sequences are in -- this was necessary once, when we actually didn't know the
    vector that had been sent down to NGS. A bit hacky.

    :param fastq_filepath:
    :param coding_startseq:
    :param noncoding_startseq:
    :param min_end_len:
    :param max_end_len:
    :param n_reads_to_sample:
    :param codingseqs_for_inference: Default None. if this is supplied, then only consider
    the read if the coding sequence matches one of codingseqs_for_inference
    exactly
    :param min_longest_support:
    :param strand:
    :return:
    """
    fastq = pysam.Fastqfile(fastq_filepath)

    start_seq_counts = dict()
    end_seq_counts = dict()

    n_sampled = 0
    while True:
        if n_sampled >= n_reads_to_sample:
            break
        try:
            aread = fastq.next()
        except Exception:
            break
        n_sampled += 1
        if codingseqs_for_inference:
            codingseq = calc_codingseq_safe(aread.sequence)
            if not codingseq or not codingseq in codingseqs_for_inference:
                continue
        if coding_startseq in aread.sequence:
            end_before_codingstart = aread.sequence[0:aread.sequence.index(coding_startseq)]
            if len(end_before_codingstart) >= min_end_len and len(end_before_codingstart) <= max_end_len:
                if end_before_codingstart not in start_seq_counts:
                    start_seq_counts[end_before_codingstart] = 0
                start_seq_counts[end_before_codingstart] += 1
        if noncoding_startseq in aread.sequence:
            end_after_noncodingstart = aread.sequence[
                                       aread.sequence.index(noncoding_startseq) + len(noncoding_startseq):]
            if len(end_after_noncodingstart) >= min_end_len and len(end_after_noncodingstart) <= max_end_len:
                if end_after_noncodingstart not in end_seq_counts:
                    end_seq_counts[end_after_noncodingstart] = 0
                end_seq_counts[end_after_noncodingstart] += 1
    most_common_startseq = ""
    commonest_count = 0
    longest_start_with_support = ''
    for seq in start_seq_counts:
        if start_seq_counts[seq] > commonest_count:
            most_common_startseq = seq
            commonest_count = start_seq_counts[seq]
        if start_seq_counts[seq] >= min_longest_support and len(seq) > len(longest_start_with_support):
            longest_start_with_support = seq
    print("start commonest count: %d" % commonest_count)
    most_common_endseq = ""
    commonest_count = 0
    longest_end_with_support = ''
    for seq in end_seq_counts:
        if end_seq_counts[seq] > commonest_count:
            most_common_endseq = seq
            commonest_count = end_seq_counts[seq]
        if end_seq_counts[seq] >= min_longest_support and len(seq) > len(longest_end_with_support):
            longest_end_with_support = seq
    print("end commonest count: %d" % commonest_count)
    print("Longest with support over %d: %s, %s" % (
    min_longest_support, longest_start_with_support, longest_end_with_support))
    return (most_common_startseq, most_common_endseq)


def fasta2gff(fasta_file, out_file):
    """Create a GFF file with a trivial GFF entry for each FASTA
       sequence, assuming that each FASTA sequence has been encoded in
       a reference sequence database under its own name.
       "fasta_proteins" is a misnomer that would need a refactor to fix
       (see module protein), since these are DNA sequences, not protein."""
    fasta_proteins = fasta_io.load_fasta_proteins(fasta_file)
    gff_source = basename(fasta_file.name)
    for fasta_protein in fasta_proteins:
        gff_seqname = fasta_protein.name
        gff_start = 1
        gff_end = len(fasta_protein.sequence)
        gff_line = make_gff_line(gff_seqname, gff_source, gff_seqname,
                                 gff_start, gff_end)
        out_file.write(gff_line + '\n')


def make_gff_line(chromosome, source, feature_name, start_1based, end_1based,
                  score='.', strand='.', frame='.', attribute=''):
    """Create a line of a GFF file.
    Other people have written plenty of ways to do this. This is quick and
    dirty. Consider switching to BCBio if we need more functionality"""
    return '\t'.join([chromosome, source, feature_name,
                      str(start_1based), str(end_1based),
                      str(score), strand, str(frame), attribute])


def write_analysis_to_file(library_analysis, outfile,
                           should_sort=True):
    """write all ObservedCodingSeqs from a sample analysis to a file, in tsv format"""
    if should_sort:
        library_analysis.observed_codingseqs.sort(key=lambda x: x.readcount, reverse=True)
    outfile.write("\t".join(['sequence', 'readcount', 'readcount_posstrand', 'has_refseq',
                             'refseq_name', 'refseq_sequence', 'n_diffs', 'diff_positions']) + '\n')
    for obscodingseq in library_analysis.observed_codingseqs:
        outfile.write(str(obscodingseq) + '\n')


def read_analysis_from_file(infile):
    """read all ObservedCodingSeqs from a library analysis file"""
    obscodingseqs = list()
    csvreader = csv.reader(infile, delimiter='\t')
    for row in csvreader:
        # ignore first row (header) if it's there
        if row[0] == 'sequence':
            continue
        sequence = row[0]
        readcount = int(row[1])
        readcount_posstrand = int(row[2])
        has_refseq = row[3]
        refseq_name = row[4]
        refseq_sequence = row[5]
        refseq_diffpositions_str = row[6]
        refseq_diffpositions = []
        if len(refseq_diffpositions_str) > 0:
            [int(x) for x in refseq_diffpositions_str.split(',')]
        obscodingseq = ObservedCodingSeq(sequence, readcount, readcount_posstrand,
                                         has_refseq == 'Y', refseq_name, refseq_sequence, refseq_diffpositions)
        obscodingseqs.append(obscodingseq)
    return HTSLibraryAnalysis(obscodingseqs)


class ObservedCodingSeq:
    """a class to store information about a sequence observed in HTS data in
    one or more reads, and to summarize the information """

    def __init__(self, sequence, readcount,
                 readcount_posstrand, has_refseq,
                 refseq_name, refseq_sequence, refseq_diffpositions):
        self.sequence = sequence
        self.readcount = readcount
        self.readcount_posstrand = readcount_posstrand
        self.has_refseq = has_refseq
        self.refseq_name = refseq_name
        self.refseq_sequence = refseq_sequence

        ###BROKEN!!!
        self.refseq_diffpositions = refseq_diffpositions
        self.mut_positions_refbased = list()
        self.ins_positions_after_refbased = list()
        self.del_positions_refbased = list()
        self.all_bad_positions_refbased = list()


    def __str__(self):
        diffpos_str = ','.join([str(x) for x in self.refseq_diffpositions])
        has_refseq_str = 'N'
        if self.has_refseq:
            has_refseq_str = 'Y'
        return '\t'.join([self.sequence, str(self.readcount),
                          str(self.readcount_posstrand),
                          has_refseq_str, self.refseq_name,
                          self.refseq_sequence, str(len(self.refseq_diffpositions)), diffpos_str])


class HTSLibraryAnalysis:
    """A class to store the results of a time-consuming analysis of a HTS
    run of a library sample. """

    def __init__(self, observed_codingseqs,
                 min_reads=1,
                 should_sort=True):
        self.observed_codingseqs = observed_codingseqs
        self.seq_obscodingseq_map = dict()
        for obscodingseq in observed_codingseqs:
            if obscodingseq.readcount >= min_reads:
                self.seq_obscodingseq_map[obscodingseq.sequence] = obscodingseq

    def get_sequences(self):
        return self.seq_obscodingseq_map.keys()

    def get_obscodingseqs(self):
        return self.seq_obscodingseq_map.values()

    # chart-building tools. A lot of charts are built multiple times
    # for different read depths

    def build_readpos_hist(self, min_readcount, max_readcount, title):
        diffpositions = list()
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >= min_readcount and obscodingseq.readcount <= max_readcount:
                diffpositions.extend(obscodingseq.all_bad_positions_refbased)
        if len(diffpositions) > 2:
            return charts.hist(diffpositions, title)
        return None

    def build_readcount_hist(self, min_readcount, max_readcount, logmode=False):
        readcounts = list()
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >= min_readcount and obscodingseq.readcount <= max_readcount:
                to_append = obscodingseq.readcount
                if logmode:
                    to_append = math.log(to_append)
                readcounts.append(to_append)
        if len(readcounts) > 2:
            title = "Read counts, " + str(min_readcount) + " to " + str(max_readcount)
            if logmode:
                title = "Log " + title
            return charts.hist(readcounts, title)
        return None

    def build_proteinlength_hist(self, min_readcount, max_readcount):
        readcounts = list()
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >= min_readcount and obscodingseq.readcount <= max_readcount:
                readcounts.append(len(dna.forward_translate_dna_firstframe(obscodingseq.sequence, stop_before_stop=True)))
        if len(readcounts) > 2:
            return charts.hist(readcounts,
                               "Coded protein lengths, " + str(min_readcount) + " to " + str(max_readcount) + " reads")
        return None

    def build_protein_list(self, min_reads, fasta_proteins,
                           min_length=MIN_OUTPUT_PROTEIN_LENGTH):
        """ build a list of proteins based on this analysis, with min_reads
        reads or more supporting at least one coding sequence for the protein """
        result = list()
        protein_seqs = set()

        fasta_proteinseq_protein_map = dict()
        for protein in fasta_proteins:
            fasta_proteinseq_protein_map[protein.sequence] = protein

        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount < min_reads:
                continue
            proteinseq = dna.forward_translate_dna_firstframe(obscodingseq.sequence, True)
            if proteinseq in protein_seqs or len(proteinseq) < min_length:
                continue
            protein_seqs.add(proteinseq)

            if proteinseq in fasta_proteinseq_protein_map:
                result.append(fasta_proteinseq_protein_map[proteinseq])
                continue

            protein_name = 'unknown_' + proteinseq
            if obscodingseq.has_refseq:
                protein_name = obscodingseq.refseq_name
                if proteinseq != dna.forward_translate_dna_firstframe(obscodingseq.refseq_sequence):
                    protein_name = "corrupted_" + protein_name + "_" + proteinseq

            result.append(proteins.Protein(protein_name, proteinseq))

        return result

    def build_proteinseq_readcount_map(self):
        """build a map from protein sequences to read counts"""
        result = dict()
        for obscodingseq in self.observed_codingseqs:
            proteinseq = dna.forward_translate_dna_firstframe(obscodingseq.sequence)
            if proteinseq in result:
                result[proteinseq] = max(obscodingseq.readcount, result[proteinseq])
            else:
                result[proteinseq] = obscodingseq.readcount
        return result

    def compare_readcounts_proteinlist_present_absent(self, obsproteinlist, allproteinlist):
        readcounts_present = []
        readcounts_absent = []
        allproteinlist_sequences = [protein.sequence for protein in allproteinlist]
        obsproteinlist_sequences = [protein.sequence for protein in obsproteinlist]

        proteinseq_readcount_map = self.build_proteinseq_readcount_map()
        for proteinseq in proteinseq_readcount_map:
            if proteinseq in allproteinlist_sequences:
                if proteinseq in obsproteinlist_sequences:
                    readcounts_present.append(proteinseq_readcount_map[proteinseq])
                else:
                    readcounts_absent.append(proteinseq_readcount_map[proteinseq])
        charts.multihist([readcounts_present, readcounts_absent],
                         title='read counts for %d proteins in list (blue), %d not in list (red)' %
                               (len(readcounts_present), len(readcounts_absent))).show()
        charts.hist(readcounts_present,
                    title='read counts for %d proteins in list (blue)' %
                          len(readcounts_present)).show()
        charts.hist(readcounts_absent,
                    title='read counts for %d proteins not in list (blue)' %
                          len(readcounts_absent)).show()

    def build_dna_fasta_entry_list(self, min_reads, base_dnaseqs):
        """build a fasta file's worth of entries for the DNA sequences
        we found."""
        fasta_entries = list()
        for obscodingseq in self.observed_codingseqs:
            name = 'unknown_' + obscodingseq.sequence
            if obscodingseq.has_refseq:
                name = obscodingseq.refseq_name
                if not obscodingseq.sequence in base_dnaseqs:
                    name = "corrupted_" + name + "_" + obscodingseq.sequence

            fasta_entries.append(proteins.Protein(name, obscodingseq.sequence))

        return fasta_entries

    def predict_proteinseqs(self, min_readcount=40, min_proteinlength=20,
                            max_proteinlength=65):
        """for casual use"""
        result = set()
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >= min_readcount:
                proteinseq = dna.forward_translate_dna_firstframe(obscodingseq.sequence, True)
                if len(proteinseq) >= min_proteinlength and len(proteinseq) <= max_proteinlength:
                    result.add(proteinseq)
        return result

    def build_hts_analysis_charts(self, pdf_file, known_sequences=None,
                                  max_unknownseqs=200,
                                  threeprime_codingstart='GGATCC',
                                  fiveprime_noncodingstart='TAATGCGGCCGC'):
        """
        Build charts summarizing this analysis. This is very hokey. Most of these are histograms, and, as you can
        see, I'm segmenting the histograms to provide more detailed information about some low-count regions
        that are hidden in the context of the whole sample. It would probably be much better to histogram the
        log space.
        :param pdf_file:
        :param known_sequences:
        :param max_unknownseqs:
        :param threeprime_codingstart:
        :param fiveprime_noncodingstart:
        :return:
        """
        pdf_charts = list()

        if known_sequences:
            known_codingseqs = [calc_codingseq_safe(seq, threeprime_codingstart, fiveprime_noncodingstart) for seq in
                                known_sequences]

        pdf_charts.append(self.build_readcount_hist(1, 999999999999999, logmode=True))

#        for readcount_range in [(1, 9999999999)]:
#            mychart = self.build_readcount_hist(readcount_range[0], readcount_range[1])
#            if mychart:
#                pdf_charts.append(mychart)

        for readcount_range in [(1, 9999999999),
                                (30, 9999999999)]:
            mychart = self.build_proteinlength_hist(readcount_range[0], readcount_range[1])
            if mychart:
                pdf_charts.append(mychart)

        n_inframe = 0
        n_outofframe = 1
        n_mapped = 0
        n_unmapped = 0
        min_readcount_for_framecheck = 3
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >=min_readcount_for_framecheck:
                if len(obscodingseq.sequence) % 3 == 0:
                    n_inframe += 1
                else:
                    n_outofframe += 1
                if obscodingseq.has_refseq:
                    n_mapped += 1
                else:
                    n_unmapped += 1
        pdf_charts.append(charts.pie([n_inframe, n_outofframe],
                                     labels=['in frame', 'out of frame'],
                                     title="In (%d) and out (%d) of frame (>= %d reads)" %
                                           (n_inframe, n_outofframe, min_readcount_for_framecheck)))
        print("In frame: %d. Out of frame: %d" % (n_inframe, n_outofframe))

        readpos_hist = self.build_readpos_hist(3, 999999999,
                                               "seq positions different (3+ reads)")
        if readpos_hist:
            pdf_charts.append(readpos_hist)

        nt_list = ['A', 'T','G','C']

        if known_sequences:
            nt_count_map = dna.calc_nt_proportion_map(known_sequences)
            pdf_charts.append(charts.pie([nt_count_map[x] for x in nt_list], labels=nt_list,
                                         title="library sequence nt proportions"))
        # show distribution of nucleotides before mutations
        nt_badcount_map = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for obscodingseq in self.observed_codingseqs:
            for badpos in obscodingseq.mut_positions_refbased:
                if badpos-1 in xrange(0, len(obscodingseq.refseq_sequence)):
                    nt_badcount_map[obscodingseq.refseq_sequence[badpos-1]] += 1
        pdf_charts.append(charts.pie([nt_badcount_map[x] for x in nt_list], labels=nt_list,
                          title="Residues before %d mutations" % sum(nt_badcount_map.values())))

        # show distribution of nucleotides before insertions
        nt_badcount_map = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for obscodingseq in self.observed_codingseqs:
            for badpos in obscodingseq.ins_positions_after_refbased:
                if badpos in xrange(0, len(obscodingseq.refseq_sequence)):
                    nt_badcount_map[obscodingseq.refseq_sequence[badpos]] += 1
        pdf_charts.append(charts.pie([nt_badcount_map[x] for x in nt_list], labels=nt_list,
                                     title="Residues before %d insertions" % sum(nt_badcount_map.values())))

        # show distribution of nucleotides before deletions
        nt_badcount_map = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for obscodingseq in self.observed_codingseqs:
            for badpos in obscodingseq.ins_positions_after_refbased:
                if badpos-1 in xrange(0, len(obscodingseq.refseq_sequence)):
                    nt_badcount_map[obscodingseq.refseq_sequence[badpos-1]] += 1
        pdf_charts.append(charts.pie([nt_badcount_map[x] for x in nt_list], labels=nt_list,
                                     title="Residues before %d deletions" % sum(nt_badcount_map.values())))

        # barchart indels vs. mutations
        n_with_indels = 0
        n_with_mutations = 0
        for obscodingseq in self.observed_codingseqs:
            if len(obscodingseq.ins_positions_after_refbased) + len(obscodingseq.del_positions_refbased) > 0:
                n_with_indels += 1
            if len(obscodingseq.mut_positions_refbased) > 0:
                n_with_mutations += 1
        pdf_charts.append(charts.multibar([[n_with_indels, n_with_mutations]], labels=["indels","mutations"]))

        ns_diffs = list()
        ns_readcounts = list()
        min_readcount = 3
        for obscodingseq in self.observed_codingseqs:
            if obscodingseq.readcount >= min_readcount:
                ns_diffs.append(len(obscodingseq.all_bad_positions_refbased))
                ns_readcounts.append(obscodingseq.readcount)
        if len(ns_diffs) > 2:
            pdf_charts.append(charts.hist(ns_diffs, "numbers of differences (3+ reads"))
            pdf_charts.append(charts.scatterplot(ns_diffs, ns_readcounts, "#diffs (x) vs. readcount (y, 3+ reads)"))

        if known_sequences:
            print(" Building known-sequence charts...")
            known_sequence_readcount_countsmap = dict()
            total_readcount_countsmap = dict()
            for obscodingseq in self.observed_codingseqs:
                if not obscodingseq.readcount in total_readcount_countsmap:
                    total_readcount_countsmap[obscodingseq.readcount] = 0
                total_readcount_countsmap[obscodingseq.readcount] += 1
                if obscodingseq.sequence in known_codingseqs:
                    if not obscodingseq.readcount in known_sequence_readcount_countsmap:
                        known_sequence_readcount_countsmap[obscodingseq.readcount] = 0
                    known_sequence_readcount_countsmap[obscodingseq.readcount] += 1

            knownseq_count_dataset = list()
            unknownseq_count_dataset = list()
            knownseq_proportion_dataset = list()
            knownseq_cumcount = 0
            allseq_cumcount = 0
            for i in sorted(total_readcount_countsmap, reverse=True):
                allseq_cumcount += total_readcount_countsmap[i]
                if i in known_sequence_readcount_countsmap:
                    knownseq_cumcount += known_sequence_readcount_countsmap[i]

                knownseq_count_dataset.append((i, knownseq_cumcount))
                knownseq_proportion_dataset.append((i, float(knownseq_cumcount) / float(allseq_cumcount)))
                unknownseq_count_dataset.append((i, min(allseq_cumcount - knownseq_cumcount, max_unknownseqs)))
            knownseq_count_dataset.reverse()
            knownseq_proportion_dataset.reverse()
            unknownseq_count_dataset.reverse()
            pdf_charts.append(charts.line_plot([mypair[0] for mypair in knownseq_count_dataset],
                                               [mypair[1] for mypair in knownseq_count_dataset],
                                               "Counts of known sequences at read counts (desc)"))
            pdf_charts.append(charts.line_plot([mypair[0] for mypair in knownseq_proportion_dataset],
                                               [mypair[1] for mypair in knownseq_proportion_dataset],
                                               "Proportions of sequences known at read counts (desc)"))
            pdf_charts.append(charts.line_plot([mypair[0] for mypair in unknownseq_count_dataset],
                                               [mypair[1] for mypair in unknownseq_count_dataset],
                                               "Counts of unknown sequences at read counts (desc, max=%d)" % max_unknownseqs))

        print("Writing pdf file %s with %d charts." % (pdf_file.name, len(pdf_charts)))
        charts.write_pdf(pdf_charts, pdf_file)
        pdf_file.close()


def build_alignment_script(fastq_file, dna_fasta_file,
                           alignment_dir, out_file,
                           read_length=DEFAULT_READ_LENGTH,
                           n_cores=DEFAULT_N_CORES,
                           email='damonmay@uw.edu',
                           max_diffs=DEFAULT_MAX_ALIGN_DIFFS,
                           library_file=None,
                           min_reads_for_output=20,
                           run_cutadapt=False,
                           build_charts=False):
    """
    Build a script to perform an alignment of an HTS run

    :param fastq_file:
    :param dna_fasta_file:
    :param alignment_dir:
    :param out_file:
    :param read_length:
    :param n_cores:
    :param email:
    :param max_diffs:
    :param library_file:
    :param min_reads_for_output:
    :param run_cutadapt:
    :param build_charts: This is actually a boolean indicator of whether to run analyze_library_hts. If False, we
    still generate the command to run, but we comment it out so that it can be run later. This is done so that
    alignment can be separated from analysis, since both of those things can take a lot of time and analysis might
    need to be re-run.
    :return:
    """
    untrimmed_fastq_basename = basename(fastq_file.name)
    if run_cutadapt:
        fastq_basename = untrimmed_fastq_basename[:untrimmed_fastq_basename.index('fastq.gz')] + 'trimmed.fastq.gz'
    else:
        fastq_basename = untrimmed_fastq_basename
    fastq_filename = fastq_basename
    if '.fastq' in fastq_basename:
        fastq_basename = fastq_basename[:fastq_basename.index('.fastq')]
    sai_filename = fastq_basename + '.sai'
    bam_filename = fastq_basename + '.bam'
    sorted_bam_basename = fastq_basename + '.sorted'
    sorted_bam_filename = fastq_basename + '.sorted.bam'

    script_lines = list()
    script_lines.append("#!/bin/bash")
    script_lines.append("#SBATCH -n%d -N1 -t 0-4 -p campus --mail-type=ALL --mail-user=%s" % (n_cores, email))
    script_lines.append("")
    script_lines.append("export PATH=/home/solexa/apps/bwa/bwa-0.7.5a/:$PATH")
    if run_cutadapt:
        script_lines.append("runcutadapt %s" % untrimmed_fastq_basename)
    script_lines.append("mkdir -p %s" % alignment_dir)
    script_lines.append("cd %s" % alignment_dir)
    script_lines.append("bwa aln -t %d -k %d -l %d %s ../%s > %s" %
                        (n_cores, max_diffs, read_length, dna_fasta_file.name,
                         fastq_filename, sai_filename))
    script_lines.append("bwa samse %s %s ../%s | samtools view -bS - > %s" %
                        (dna_fasta_file.name, sai_filename, fastq_filename,
                         bam_filename))
    script_lines.append("samtools flagstat %s" % bam_filename)
    script_lines.append("rm %s" % sai_filename)

    script_lines.append("echo sorting")
    script_lines.append("samtools sort %s %s" % (bam_filename, sorted_bam_basename))
    script_lines.append("samtools index %s" % (sorted_bam_filename))

    out_file.write("\n".join(script_lines))
    out_file.write("\n")
    out_file.flush()
    out_file.close()


def assess_quals(samfile, n_reads_to_eval=500000):
    """histogram the quality scores of mapped and unmapped positions"""
    quals_mappedbases = list()
    quals_unmappedbases = list()
    n_evaled = 0
    for aread in samfile:
        if n_evaled >= n_reads_to_eval:
            break
        n_evaled += 1
        quals = [ord(x) for x in aread.qual]
        for i in xrange(0, len(aread.positions)):
            if aread.positions[i]:
                quals_mappedbases.append(quals[i])
            else:
                quals_unmappedbases.append(quals[i])
    print("Quals: unmapped: %d, mapped: %d" % (len(quals_unmappedbases), len(quals_mappedbases)))
    charts.hist(quals_mappedbases).show()

    print()

        
    
        
    
