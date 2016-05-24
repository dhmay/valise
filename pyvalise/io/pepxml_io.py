#!/usr/bin/env python

"""
PepXML reading and (eventually) writing
# todo: probably does not handle terminal modifications on IDs correctly.
"""

import logging
import re

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
from pyvalise.proteomics import peptides
import xml.dom.minidom as minidom


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2012-2014 Fred Hutchinson Cancer Research Center"
__license__ = ""
__version__ = ""

PEPXML_NS_URL = "http://regis-web.systemsbiology.net/pepXML"
PEPXML_NS = "{" + PEPXML_NS_URL + "}"

# minimum ratio, for pegging 0-ratio values
MIN_RATIO = 0.001

log = logging.getLogger(__name__)


def read_pepxml(pepxml_file,
                min_pprophet=None):
    """Converts a pepXML file into a MSMSPipelineAnalysis.
       Only includes a PeptideIdentification if PeptideProphet probability
       >= min_pprophet_probability (if specified
    """

    #can't use pyteomics for this. Doesn't handle my pepxml files for some reason
    tree = ET.parse(pepxml_file)
    root = tree.getroot()
    url = re.search(r"\{.*\}", root.tag).group(0)

    if min_pprophet is not None:
        log.debug("Peptide probability cutoff: %f" % min_pprophet)

    result = peptides.MSMSPipelineAnalysis()

    for runsummary_elem in root.findall(url + 'msms_run_summary'):
        log.debug("read_pepxml found run.")

        searchsummary_elem = runsummary_elem.find(url + 'search_summary')
        search_engine = searchsummary_elem.get('search_engine')
        fasta = searchsummary_elem.find(url + 'search_database').get('local_path')
        max_missed_cleavages = searchsummary_elem.find(url + 'enzymatic_search_constraint').get(
            'max_num_internal_cleavages')
        min_tryptic_termini = searchsummary_elem.find(url + 'enzymatic_search_constraint').get('min_number_termini')
        run_mods = []

        for aamod_elem in searchsummary_elem.findall(url + 'aminoacid_modification'):
            run_mods.append(peptides.AminoacidModification(aamod_elem.get('aminoacid'),
                                                           float(aamod_elem.get('massdiff')),
                                                           aamod_elem.get('variable') == 'Y'))
        for term_modmod_elem in searchsummary_elem.findall(url + 'terminal_modification'):
            run_mods.append(peptides.AminoacidModification(term_modmod_elem.get('terminus'),
                                                           float(term_modmod_elem.get('massdiff')),
                                                           term_modmod_elem.get('variable') == 'Y'))
        run_result = peptides.MSRunSearchResult(
            search_engine=search_engine,
            max_missed_cleavages=max_missed_cleavages,
            min_tryptic_termini=min_tryptic_termini,
            modifications=run_mods,
            name=runsummary_elem.attrib["base_name"],
            fasta=fasta)

        result.runs.append(run_result)

        for spectrum_query_elem in runsummary_elem.findall(url + 'spectrum_query'):
            scan = int(spectrum_query_elem.attrib["start_scan"])
            time = float(spectrum_query_elem.attrib["retention_time_sec"])
            charge = int(spectrum_query_elem.attrib["assumed_charge"])
            observed_mass = float(spectrum_query_elem.get('precursor_neutral_mass'))
            spectrum_name = spectrum_query_elem.get('spectrum')

            # assuming only one search_result per spectrum_query and
            # only one search_hit per search_result (or we only care about
            # the first)
            search_hit_elem = spectrum_query_elem.find(url + 'search_result').find(url + 'search_hit')
            if not search_hit_elem:
                continue
            peptide = search_hit_elem.get("peptide")
            prev_aa = search_hit_elem.get('peptide_prev_aa')
            next_aa = search_hit_elem.get('peptide_next_aa')
            num_tol_term = int(search_hit_elem.get('num_tol_term'))
            mass = float(search_hit_elem.get("calc_neutral_pep_mass"))
            proteins = list()
            proteins.append(search_hit_elem.attrib["protein"])
            for alt_protein_elem in search_hit_elem.findall(url + 'alternative_protein'):
                proteins.append(alt_protein_elem.attrib["protein"])
            probability = 0
            for analysis_result_elem in search_hit_elem.findall(url + 'analysis_result'):
                # This is hacky. Shouldn't be a loop
                for peptideprophet_result_elem in analysis_result_elem.findall(url + 'peptideprophet_result'):
                    probability = float(peptideprophet_result_elem.attrib["probability"])
            modifications = list()
            modinfo_elem = search_hit_elem.find(url + 'modification_info')
            if modinfo_elem:
                for mod_elem in modinfo_elem.findall(url + 'mod_aminoacid_mass'):
                    # positions in pepXML are 1-based
                    position = int(mod_elem.get('position')) - 1
                    aa_at_position = peptide[position]
                    modified_mass = float(mod_elem.get('mass'))
                    modifications.append(peptides.ModifiedAminoacid(
                        position,
                        modified_mass,
                        modified_mass - peptides.get_aa_mass(aa_at_position)))

            # read the heavy/light ratio
            ratio_heavy_light = None
            quant_heavy_area = None
            quant_light_area = None
            quant_labelfree_peakintensity = None
            quant_labelfree_peakarea = None
            quant_labelfree_peak_rt_seconds = None
            quant_labelfree_start_rt_seconds = None
            quant_labelfree_end_rt_seconds = None
            for analysis_elem in search_hit_elem.findall(url + 'analysis_result'):
                # only XPress quantitation supported now
                if analysis_elem.attrib['analysis'] == 'xpresslabelfree':
                    xpress_elem = analysis_elem.find(url + 'xpresslabelfree_result')
                    if xpress_elem is not None:
                        quant_labelfree_peakintensity = float(xpress_elem.attrib['peak_intensity'])
                        quant_labelfree_peakarea = float(xpress_elem.attrib['peak_area'])
                        quant_labelfree_peak_rt_seconds = float(xpress_elem.attrib['peak_intensity_RT_seconds'])
                        quant_labelfree_start_rt_seconds = float(xpress_elem.attrib['first_scan_RT_seconds'])
                        quant_labelfree_end_rt_seconds = float(xpress_elem.attrib['last_scan_RT_seconds'])
                    break
                elif analysis_elem.attrib['analysis'] == 'xpress':
                    xpress_elem = analysis_elem.find(url + 'xpressratio_result')
                    if xpress_elem is not None:
                        ratio_heavy_light = 1.0 / max(float(xpress_elem.attrib['decimal_ratio']), MIN_RATIO)
                        quant_heavy_area = float(xpress_elem.attrib['heavy_area'])
                        quant_light_area = float(xpress_elem.attrib['light_area'])
                    break

            if not min_pprophet or probability >= min_pprophet:
                peptide_id = peptides.PeptideIdentification(scan, time,
                                                            peptide, probability,
                                                            mass, charge,
                                                            proteins,
                                                            observed_mass=observed_mass,
                                                            prev_aa=prev_aa,
                                                            next_aa=next_aa,
                                                            modifications=modifications,
                                                            spectrum_name=spectrum_name,
                                                            num_tol_term=num_tol_term,
                                                            ratio_heavy_light=ratio_heavy_light,
                                                            quant_heavy_area=quant_heavy_area,
                                                            quant_light_area=quant_light_area,
                                                            quant_labelfree_peakintensity=quant_labelfree_peakintensity,
                                                            quant_labelfree_peakarea=quant_labelfree_peakarea,
                                                            quant_labelfree_peak_rt_seconds=quant_labelfree_peak_rt_seconds,
                                                            quant_labelfree_start_rt_seconds=quant_labelfree_start_rt_seconds,
                                                            quant_labelfree_end_rt_seconds=quant_labelfree_end_rt_seconds)
                run_result.peptide_ids.append(peptide_id)

    return result


def write_pepxml(msms_pipeline_analysis, outfile):
    """
    takes a peptides.MSMSPipelineAnalysis and writes it out to a file.
    """
    # TODO: write ratios
    ET.register_namespace('', PEPXML_NS_URL)
    root = ET.Element(PEPXML_NS + "msms_pipeline_analysis")


    #add dummy PeptideProphet info
    ET.SubElement(
        ET.SubElement(root, PEPXML_NS + 'analysis_summary',
                      {'analysis': 'peptideprophet'}),
        PEPXML_NS + 'peptideprophet_summary', {'version': 'PeptideProphet (unknown)',
                                               'min_prob': '0.05'})

    for run in msms_pipeline_analysis.runs:
        run_elem = ET.SubElement(root, PEPXML_NS + "msms_run_summary",
                                 {'base_name': run.name})
        # add sample enzyme info
        # todo: make generic beyond trypsin
        add_trypsin_info(run_elem)

        search_summary_elem = ET.SubElement(run_elem, PEPXML_NS + "search_summary",
                                            {'precursor_mass_type': 'monoisotopic',
                                             'search_engine': run.search_engine,
                                             'base_name': run.name})
        fasta = run.fasta
        if fasta is None:
            fasta = ''
        ET.SubElement(search_summary_elem, PEPXML_NS + 'search_database',
                      {'local_path': fasta})
        ET.SubElement(search_summary_elem, PEPXML_NS + 'enzymatic_search_constraint',
                      {'max_num_internal_cleavages': str(run.max_missed_cleavages),
                       'min_number_termini': str(run.min_tryptic_termini),
                       'enzyme': 'trypsin'})
        for mod in run.modifications:
            if mod.is_terminal():
                # todo: figure out what I ougtha be doing with protein_terminus
                ET.SubElement(search_summary_elem, PEPXML_NS + 'terminal_modification',
                              {'mass': str(mod.modified_mass),
                               'massdiff': str(mod.modified_mass),
                               'variable': mod.get_variable_Y_N(),
                               'terminus': mod.aa,
                               'protein_terminus': 'N'})
            else:
                # todo: handle symbol "properly" if necessary
                ET.SubElement(search_summary_elem, PEPXML_NS + 'aminoacid_modification',
                              {'aminoacid': mod.aa,
                               'massdiff': str(mod.massdiff),
                               'mass': str(mod.modified_mass),
                               'variable': mod.get_variable_Y_N(),
                               'symbol': '^'})
        #print([run.name, run.search_engine, run.fasta])
        #print(ET.tostring(run_elem))
        sq_index = 0
        for peptide_id in run.peptide_ids:
            sq_index += 1
            add_spectrum_query(run_elem, peptide_id, sq_index, run.name)

    xml = minidom.parseString(ET.tostring(root))
    outfile.write(xml.toprettyxml(indent='  '))
    outfile.close()
    log.debug("write.pepxml: done")


def add_spectrum_query(run_elem, peptide_id, index, run_base_name):
    """ add a spectrum_query to an msms_run_summary.
          index: 1-based number for each spectrum_query. """
    prec_neut_mass = peptide_id.mass
    if peptide_id.observed_mass:
        prec_neut_mass = peptide_id.observed_mass
    else:
        peptide_id.observed_mass = peptide_id.mass
    if not peptide_id.num_tol_term:
        peptide_id.num_tol_term = 2
    # print([index, peptide_id.scan, peptide_id.charge, run_base_name, peptide_id, prec_neut_mass, peptide_id.time])
    sq_elem = ET.SubElement(run_elem, PEPXML_NS + 'spectrum_query',
                            {'index': str(index),
                             'start_scan': str(peptide_id.scan),
                             'end_scan': str(peptide_id.scan),
                             'assumed_charge': str(peptide_id.charge),
                             'spectrum': construct_spectrum_name(run_base_name, peptide_id),
                             'precursor_neutral_mass': str(prec_neut_mass),
                             'retention_time_sec': str(peptide_id.time)})

    prev_aa_string = '-'
    if peptide_id.prev_aa:
        prev_aa_string = peptide_id.prev_aa
    next_aa_string = '-'
    if peptide_id.next_aa:
        next_aa_string = peptide_id.next_aa
    protein = ""
    if peptide_id.proteins:
        protein = peptide_id.proteins[0]
    # todo: num_matched_ions is fake
    # todo: num_missed_cleavages is fake
    # print([peptide_id.sequence, peptide_id.observed_mass, len(peptide_id.proteins), peptide_id.get_delta_mass(),
    #        peptide_id.calc_n_missed_cleavages(), peptide_id.num_tol_term, prev_aa_string, next_aa_string])

    searchhit_elem = ET.SubElement(ET.SubElement(sq_elem, PEPXML_NS + 'search_result'),
                                   PEPXML_NS + 'search_hit',
                                   {'peptide': peptide_id.sequence,
                                    'calc_neutral_pep_mass': str(peptide_id.observed_mass),
                                    'hit_rank': '1',
                                    'num_matched_ions': '5',
                                    'num_tot_proteins': str(len(peptide_id.proteins)),
                                    'massdiff': str(peptide_id.get_delta_mass()),
                                    'num_missed_cleavages': str(peptide_id.calc_n_missed_cleavages()),
                                    'is_rejected': '0',
                                    'num_tol_term': str(peptide_id.num_tol_term),
                                    'peptide_prev_aa': prev_aa_string,
                                    'peptide_next_aa': next_aa_string,
                                    'protein': protein})
    if peptide_id.proteins and len(peptide_id.proteins) > 1:
        for protein in peptide_id.proteins[1:]:
            ET.SubElement(searchhit_elem, PEPXML_NS + 'alternative_protein',
                          {'protein': protein,
                           'protein_descr': 'DUMMY'})
    if peptide_id.modifications:
        modinfo_elem = ET.SubElement(searchhit_elem, PEPXML_NS + 'modification_info')
        for mod in peptide_id.modifications:
            ET.SubElement(modinfo_elem, PEPXML_NS + 'mod_aminoacid_mass',
                          {'position': str(mod.get_position_1based()),  # positions in pepXML are 1-based
                           'mass': str(mod.modified_mass)})
    # todo: add search scores
    # todo: figure out all_ntt_prob
    #print(peptide_id.probability)
    ET.SubElement(ET.SubElement(searchhit_elem, PEPXML_NS + 'analysis_result',
                                {'analysis': 'peptideprophet'}),
                  PEPXML_NS + 'peptideprophet_result',
                  {'probability': str(peptide_id.probability),
                   'all_ntt_prob': '(0.0000,0.0000,%f)' % peptide_id.probability})
    #print(ET.tostring(sq_elem))


def construct_spectrum_name(run_base_name, peptide_id):
    return ".".join([run_base_name, str(peptide_id.scan),
                     str(peptide_id.scan), str(peptide_id.charge)])


def add_trypsin_info(run_elem):
    """ add trypsin enzyme info to an msms_run_summary """
    sample_enzyme_elem = ET.SubElement(run_elem, PEPXML_NS + "sample_enzyme",
                                       {'name': 'trypsin'})
    ET.SubElement(sample_enzyme_elem, PEPXML_NS + "specificity",
                  {'cut': 'KR', 'no_cut': 'P', 'sense': 'C'})
    

