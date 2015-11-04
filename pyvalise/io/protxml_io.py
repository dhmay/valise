#!/usr/bin/env python

"""
Cribbed originally from https://github.com/boscoh/pepto/
Edited heavily.

ProtXML reading and (eventually) writing
"""

import xml.etree.ElementTree as etree
import re
import logging


logger = logging.getLogger(__name__)


def parse_attrib(attrib):
    """Simple converter of XML attrib data-strings into a dict of live
       Python variables using various re matches.
    """
    result = {}
    for key, value in attrib.items():
        if key.startswith("is_") or 'flag' in key:
            result[key] = value == 'Y'
        elif re.search(r'^\d+$', value):
            result[key] = int(value)
        elif re.search(r'^\d*[.]\d*$', value):
            result[key] = float(value)
        else:
            result[key] = value
    return result


def read_protxml(prot_xml, protein_probability_cutoff=None):
    """Converts a protein.xml file into a list of protein groups.
       Only includes a group if there's at least one protein passing
       protein_probability_cutoff (if specified)
    """
    tree = etree.parse(prot_xml)
    root = tree.getroot()
    url = re.search(r"\{.*\}", root.tag).group(0)

    if protein_probability_cutoff is not None:
        logger.debug("Protein probability cutoff: %f" % protein_probability_cutoff)

    protxml_groups = list()
    for group_elem in root.findall(url + 'protein_group'):

        group_attrs = parse_attrib(group_elem.attrib)
        group = ProtXmlGroup(group_attrs['group_number'], group_attrs['probability'])
        logger.debug('Protein group %d' % group.number)
        for protein_elem in group_elem.findall(url + 'protein'):
            protein_attrs = parse_attrib(protein_elem.attrib)

            probability = protein_attrs['probability']
            percent_coverage = 0.
            if 'percent_coverage' in protein_attrs:
                percent_coverage = protein_attrs['percent_coverage']

            if protein_probability_cutoff is not None:
                if probability < protein_probability_cutoff:
                    continue

            annotation_elem = protein_elem.find(url + 'annotation')
            description = ""
            if annotation_elem is not None:
                description = annotation_elem.attrib['protein_description']

            alt_names = list()
            for alt_protein in protein_elem.findall(url + 'indistinguishable_protein'):
                alt_names.append(alt_protein.attrib['protein_name'])

            protein_name = protein_attrs['protein_name']
            logger.debug("    Protein %s passed cutoff" % protein_name)

            peptide_seqs = set()
            for peptide_elem in protein_elem.findall(url + 'peptide'):
                peptide_attrs = parse_attrib(peptide_elem.attrib)
                peptide_seqs.add(peptide_attrs['peptide_sequence'])
            spectral_count = protein_attrs['total_number_peptides']
            protein = ProtXmlProtein(protein_name, alt_names,
                                             peptide_seqs,
                                             probability,
                                             percent_coverage,
                                             spectral_count,
                                             description,
                                             protein_attrs['n_indistinguishable_proteins'])
            group.add_protein(protein)

            for analysis_elem in protein_elem.findall(url + 'analysis_result'):
                # only XPress quantitation supported now
                if analysis_elem.attrib['analysis'] == 'xpress':
                    xpress_elem = analysis_elem.find(url + 'XPressRatio')
                    if xpress_elem is not None:
                        protein.ratio_light_heavy = float(xpress_elem.attrib['ratio_mean'])

        if group.proteins:
            protxml_groups.append(group)
    return protxml_groups


class ProtXmlGroup:
    """Stores what we need from a ProtXml protein group. Add more as needed"""

    def __init__(self, number, probability):
        self.number = number
        self.proteins = list()
        self.probability = probability

    def add_protein(self, protein):
        self.proteins.append(protein)

    def __str__(self):
        return "ProtXmlGroup %d, %d proteins" % (self.number, len(self.proteins))

    def get_protein_names(self):
        """list the protein names in this group"""
        protein_names = list()
        for protein in self.proteins:
            protein_names.append(protein.name)
        return protein_names


class ProtXmlProtein:
    """Stores what we need from a ProtXml protein. Add more as needed"""

    def __init__(self, name, alt_names, peptide_seqs, probability, percent_coverage,
                 spectral_count, description, n_indistinguishable_proteins):
        self.name = name
        self.alt_names = alt_names
        self.peptide_seqs = peptide_seqs
        self.probability = probability
        self.percent_coverage = percent_coverage
        self.spectral_count = spectral_count
        self.description = description
        self.n_indistinguishable_proteins = n_indistinguishable_proteins
        self.ratio_light_heavy = None

    def __str__(self):
        return "ProtXmlProtein: %s, prob=%f, peptides=%d, count=%d" % (self.name,
                                                                       self.probability,
                                                                       len(self.peptide_seqs),
                                                                       self.spectral_count)
    
    

