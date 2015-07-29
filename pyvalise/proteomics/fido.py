#!/usr/bin/env python

"""
Utilities for dealing with FIDO output
"""

import logging
import math
import csv

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

class FidoPeptideSequence:
    """
    A representation of a peptide sequence as expressed in FIDO. Associated with one or more scans.
    Each scan has a charge
    """
    def __init__(self, sequence, scan_charge_map):
        self.sequence = sequence
        self.scan_charge_map = scan_charge_map


class FidoProtein:
    """
    A representation of a Fido-style protein result
    """
    def __init__(self, protein_id, fido_peptides, percolator_score, percolator_rank, percolator_qvalue):
        self.protein_id = protein_id
        self.fido_peptides = fido_peptides
        self.percolator_score = percolator_score
        self.percolator_rank = percolator_rank
        self.percolator_qvalue = percolator_qvalue
        self.total_scancount = 0
        for fido_peptide in fido_peptides:
            self.total_scancount += len(fido_peptide.scan_charge_map)

def parse_fido_output(fido_outputfile):
    """

    :param fido_outputfile:
    :return: a list of FidoProteins
    """
    reader = csv.DictReader(fido_outputfile, delimiter='\t')
    logger.debug("Fields in FIDO file:")
    logger.debug("%s" % reader.fieldnames)
    fido_proteins = []
    for row in reader:
        protein_id = row['protein id']
        percolator_score = float(row['percolator score'])
        percolator_rank = float(row['percolator rank'])
        percolator_qvalue = float(row['percolator q-value'])
        pepseq_scan_charge_map = {}
        for peptide_chunk in row['peptides'].split(','):
            pepseq, scan_and_charge = peptide_chunk.split('-')
            if pepseq not in pepseq_scan_charge_map:
                pepseq_scan_charge_map[pepseq] = {}
            scan, charge = [int(x) for x in scan_and_charge.split('.')]
            pepseq_scan_charge_map[pepseq][scan] = charge
        fido_peptides = []
        for pepseq in pepseq_scan_charge_map:
            fido_peptides.append(FidoPeptideSequence(pepseq, pepseq_scan_charge_map[pepseq]))
        fido_proteins.append(FidoProtein(protein_id, fido_peptides,
                                         percolator_score, percolator_rank, percolator_qvalue))
    return fido_proteins

