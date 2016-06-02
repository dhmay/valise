#!/usr/bin/env python

"""
module for getting information from GO. Adapted from Uli Koehler's code, here:
https://techoverflow.net/blog/2013/11/18/a-geneontology-obo-v1.2-parser-in-python/
"""

import logging
from collections import defaultdict
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2015 Damon May"
__license__ = ""
__version__ = ""

logger = logging.getLogger(__name__)

NAMESPACE_MOLECULAR_FUNCTION = 'molecular_function'
NAMESPACE_BIOLOGICAL_PROCESS = 'biological_process'
NAMESPACE_CELLULAR_LOCATION = 'cellular_location'


def process_go_term(go_term):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(go_term) #Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.iteritems():
        if len(value) == 1:
            ret[key] = value[0]
    parent_ids = []
    if 'is_a' in ret:
        for is_a_value in ret['is_a']:
            parent_id = is_a_value
            if ' ' in parent_id:
                parent_id = parent_id[0:parent_id.index(' ')]
            parent_ids.append(parent_id)

    return GOTerm(ret['id'], ret['name'], ret['namespace'], parent_ids)


def parse_go_obo(obo_file):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        filename: The filename to read
    """
    current_go_term = None
    for line in obo_file:
        line = line.strip()
        if not line: continue #Skip empty
        if line == "[Term]":
            if current_go_term: yield process_go_term(current_go_term)
            current_go_term = defaultdict(list)
        elif line == "[Typedef]":
            #Skip [Typedef sections]
            current_go_term = None
        else: #Not [Term]
            #Only process if we're inside a [Term] environment
            if current_go_term is None: continue
            key, sep, val = line.partition(":")
            current_go_term[key].append(val.strip())
    #Add last term
    if current_go_term is not None:
        yield process_go_term(current_go_term)


def build_goterm_dict_from_obo(obo_file):
    """
    Build a single dict of id -> GOTerm from an obo file, with all three namespaces
    :param obo_file:
    :return:
    """
    result = {}
    for goterm in parse_go_obo(obo_file):
        result[goterm.id] = goterm
    return result


def build_namespace_goterm_dicts_from_obo(obo_file, keep_molecular_function=True,
                           keep_biological_process=True,
                           keep_cellular_location=True):
    """
    Build dicts of id -> GOTerm for each namespace that's specified. Return a map from
    namespace name to each dict
    :param obo_file:
    :param keep_molecular_function:
    :param keep_biological_process:
    :param keep_cellular_location:
    :return: a map from each GO namespace to a map from GO ids to GOTerms
    """
    result = {}
    for goterm in parse_go_obo(obo_file):
        if not keep_molecular_function and goterm.namespace == NAMESPACE_MOLECULAR_FUNCTION:
            continue
        if not keep_biological_process and goterm.namespace == NAMESPACE_BIOLOGICAL_PROCESS:
            continue
        if not keep_cellular_location and goterm.namespace == NAMESPACE_CELLULAR_LOCATION:
            continue
        if goterm.namespace not in result:
            result[goterm.namespace] = {}
        result[goterm.namespace][goterm.id] = goterm
    return result


class GOTerm:
    """
    Represents a GO term
    """

    def __init__(self, id, name, namespace, parent_ids):
        self.id = id
        self.parent_ids = parent_ids
        self.name = name
        self.namespace = namespace
