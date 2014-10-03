#!/usr/bin/env python
"""utilities for commandline tools"""

import logging
import os
import ntpath
import exceptions



log = logging.getLogger(__name__)


def argtype_dir_r(prospective_dir):
    """ defines an argparse arg that's a readable directory """
    if not os.path.isdir(prospective_dir):
        raise Exception("readable_dir:{0} is not a valid directory path".format(prospective_dir))
    if os.access(prospective_dir, os.R_OK):
        return prospective_dir
    else:
        raise Exception("readable_dir:{0} is not a readable dir".format(prospective_dir))


def argtype_dir_rw(prospective_dir):
    """ defines an argparse arg that's a readable/writeable directory """
    if not os.path.isdir(prospective_dir):
        raise Exception("readable_dir:{0} is not a valid directory path".format(prospective_dir))
    if os.access(prospective_dir, os.R_OK) and os.access(prospective_dir, os.W_OK):
        return prospective_dir
    else:
        raise Exception("readable_dir:{0} is not a readable/writable dir".format(prospective_dir))

def argtype_deltamass(deltamass_str):
    """ defines an argparse arg that's a deltamass string, must end in da or ppm
    :param deltamass_str:
    :return: a map containing keys 'delta_mass' and 'is_ppm'
    """
    deltamass_str = deltamass_str.lower()
    massval = 0
    is_ppm = False
    if deltamass_str.endswith('da'):
        massval = float(deltamass_str[:deltamass_str.index('da')])
    elif deltamass_str.endswith('ppm'):
        massval = float(deltamass_str[:deltamass_str.index('ppm')])
        is_ppm = True
    else:
        raise Exception("deltamass:{0} must end with da or ppm")
    return {
        'delta_mass': massval,
        'is_ppm': is_ppm
    }

def argtype_shiftmass(shiftmass_str):
    """ defines an argparse arg that's a shiftamass string, must end in da (only Delton unit is used)
    :param shiftamass_str:
    :return: massval float value
    """
    shiftmass_str = shiftmass_str.lower()
    if shiftmass_str.endswith('da'):
        massval = float(shiftmass_str[:shiftmass_str.index('da')])
    else:
        raise Exception("shiftmass:{0} must end with da!")
    return massval




def create_file_like_file(template_file, output_dir, output_extension):
    """Utility method for creating a new file with a name like an existing file"""
    templatename = ntpath.basename(template_file.name)
    if templatename.rfind(".") >= 0:
        templatename = templatename[0:templatename.rfind(".") + 1]
    #templatename now has the full name of the template_file up to and including
    #the last '.'
    return file(os.path.join(output_dir,
                             "".join([templatename, output_extension])), 'w')


def ensure_exactly_one_arg_present(args, name1, name2):
    """ensure that arg1 XOR arg2 is non-None. If that's not the case,
    raise an exception"""
    if (vars(args)[name1] is None) == (vars(args)[name2] is None):
        raise Exception("Either %s or %s must be specified, but not both." %
                        (name1, name2))


class ArgParseException(exceptions.Exception):
    """generic exception from parsing an argument"""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
