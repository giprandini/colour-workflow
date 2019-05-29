# -*- coding: utf-8 -*-
from aiida.orm.calculation.job.quantumespresso.simple import SimpleCalculation
from aiida.orm.data.parameter import ParameterData
from aiida.parsers.parser import Parser
from aiida.parsers.plugins.quantumespresso import QEOutputParsingError
#from aiida.parsers.plugins.quantumespresso import parse_raw_out_basic
from aiida.orm.data.array.xy import XyData
import numpy as np
from aiida.common.exceptions import InvalidOperation
from aiida.common.datastructures import calc_states


__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."


class SimpleParser(Parser):
    """
    This class is the implementation of the Parser class for Simple.
    """
   
    def __init__(self, calculation):
        """
        Initialize the instance of DosParser
        """
        # check for valid input
        if not isinstance(calculation, SimpleCalculation):
            raise QEOutputParsingError("Input calc must be a DosCalculation")

        self._calc = calculation

        super(SimpleParser, self).__init__(calculation)

    def parse_with_retrieved(self, retrieved):
        """
        Parses the datafolder, stores results.
        Retrieves simple output, and some basic information from the
        out_file, such as warnings and wall_time
        """

        # suppose at the start that the job is successful
        successful = True
        new_nodes_list = []

        # check if I'm not to overwrite anything
        state = self._calc.get_state()
        if state != calc_states.PARSING:
            raise InvalidOperation("Calculation not in {} state")

        try:
            out_folder = self._calc.get_retrieved_node()
        except KeyError:
            self.logger.error("No retrieved folder found")
            return successful, new_nodes_list

        # Read standard out
        try:
            filpath = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
            with open(filpath, 'r') as fil:
                    out_file = fil.readlines()
        except OSError:
            self.logger.error("Standard output file could not be found.")
            successful = False
            return successful, new_nodes_list

        successful = False
        for i in range(len(out_file)):
            line = out_file[-i]
            if "Matrix elements computed and saved" in line:
                successful = True
                break
        if not successful:
            self.logger.error("Computation did not finish properly")
            return successful, new_nodes_list

        parsed_data = parse_raw_out_basic(out_file, "simple")
        output_params = ParameterData(dict=parsed_data)
        # Adds warnings
        for message in parsed_data['warnings']:
            self.logger.error(message)
        # Create New Nodes List
        new_nodes_list.append( ('output_parameters', output_params) ) # output_parameters
        #new_nodes_list = [(self.get_linkname_outparams(), output_params),
        #                  (self.get_linkname_dos(), xy_data)]
        return successful,new_nodes_list



def parse_raw_out_basic(out_file, calc_name):
    """
    A very simple parser for the standard out, usually aiida.out. Currently
    only parses basic warnings and the walltime.
    :param out_file: the standard out to be parsed
    :param calc_name: the name of the calculation, e.g. PROJWFC
    :return: parsed_data
    """

    # read file
    parsed_data = {}
    parsed_data['warnings'] = []
    # critical warnings: if any is found, the calculation status is FAILED
    critical_warnings = {'Maximum CPU time exceeded':'Maximum CPU time exceeded',
                         '%%%%%%%%%%%%%%':None,
                         }

    minor_warnings = {'Warning:':None,
                      'DEPRECATED:':None,
                      }
    all_warnings = dict(critical_warnings.items() + minor_warnings.items())
    for count in range (len(out_file)):
        line = out_file[count]
        # parse the global file, for informations that are written only once
        if 'TOTAL NUMBER OF OPTIMAL BASIS VECTORS :' in line:
            parsed_data['number_optimal_basis_vectors'] = int(line.split(':')[-1])
        if calc_name in line and 'WALL' in line:
            try:
                time = line.split('CPU')[1].split('WALL')[0]
                cpu_time = line.split(':')[1].split('CPU')[0]
                parsed_data['wall_time'] = time
                parsed_data['cpu_time'] = cpu_time
            except ValueError:
                parsed_data['warnings'].append('Error while parsing wall time.')
            try:
                parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
                parsed_data['cpu_time_seconds'] = convert_qe_time_to_sec(cpu_time)
            except ValueError:
                raise QEOutputParsingError("Unable to convert wall_time in seconds.")
            # Parsing of errors
        elif any( i in line for i in all_warnings):
            message = [ all_warnings[i] for i in all_warnings.keys() if i in line][0]
            if message is None:
                message = line
            if '%%%%%%%%%%%%%%' in line:
                message  = None
                messages = parse_QE_errors(out_file,count,parsed_data['warnings'])

            # if it found something, add to log
            try:
                parsed_data['warnings'].extend(messages)
            except UnboundLocalError:
                pass
            if message is not None:
                parsed_data['warnings'].append(message)

    return parsed_data

def convert_qe_time_to_sec(timestr):
    """
    Given the walltime string of Quantum Espresso, converts it in a number of
    seconds (float).
    """
    rest = timestr.strip()

    if 'd' in rest:
        days, rest = rest.split('d')
    else:
        days = '0'

    if 'h' in rest:
        hours, rest = rest.split('h')
    else:
        hours = '0'

    if 'm' in rest:
        minutes, rest = rest.split('m')
    else:
        minutes = '0'

    if 's' in rest:
        seconds, rest = rest.split('s')
    else:
        seconds = '0.'

    if rest.strip():
        raise ValueError("Something remained at the end of the string '{}': '{}'"
                         .format(timestr, rest))

    num_seconds = (
        float(seconds) + float(minutes) * 60. +
        float(hours) * 3600. + float(days) * 86400.)

    return num_seconds

def parse_QE_errors(lines,count,warnings):
    """
    Parse QE errors messages (those appearing between some lines with
    ``'%%%%%%%%'``)
    :param lines: list of strings, the output text file as read by readlines()
    or as obtained by data.split('\n') when data is the text file read by read()
    :param count: the line at which we identified some ``'%%%%%%%%'``
    :param warnings: the warnings already parsed in the file
    :return messages: a list of QE error messages
    """

    # find the indices of the lines with problems
    found_endpoint = False
    init_problem = count
    for count2,line2 in enumerate(lines[count+1:]):
        end_problem = count + count2 + 1
        if "%%%%%%%%%%%%" in line2:
            found_endpoint = True
            break
    messages = []
    if found_endpoint:
        # build a dictionary with the lines
        prob_list = lines[init_problem:end_problem+1]
        irred_list = list(set(prob_list))
        for v in prob_list:
            if ( len(v)>0 and (v in irred_list and v not in warnings) ):
                messages.append(irred_list.pop(irred_list.index(v)))

    return messages


