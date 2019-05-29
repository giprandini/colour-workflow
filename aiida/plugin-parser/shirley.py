# -*- coding: utf-8 -*-
from aiida.orm.calculation.job.shirley import ShirleyCalculation
from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError
from aiida.parsers.plugins.quantumespresso import QEOutputParsingError
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.array.xy import XyData
from aiida.common.datastructures import calc_states
from aiida.common.exceptions import InvalidOperation

__authors__ = "The AiiDA team."
__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"

class ShirleyParser(Parser):
    """
    Shirley output parser.
    """    
    
    def __init__(self,calculation):
        """
        Initialize the instance of ShirleyParser
        """
        # check for valid input
        if not isinstance(calculation,ShirleyCalculation):
            raise OutputParsingError("Input calculation must be a "
                                     "ShirleyCalculation")
        super(ShirleyParser, self).__init__(calculation)
        
        self._eps_re_array_linkname = 'array_eps_re'
        self._eps_im_array_linkname = 'array_eps_im'
            
    def parse_with_retrieved(self, retrieved):
        """
        Parses the datafolder, stores results.
        This parser for this simple code does simply store in the DB a node
        representing the file of forces in real space
        """
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
            if "Optical properties computed and saved" in line:
                successful = True
                break
        if not successful:
            self.logger.error("Computation did not finish properly")
            return successful, new_nodes_list

        # Parse standard output
        parsed_data = parse_raw_out_basic(out_file, "simple_ip")
        output_params = ParameterData(dict=parsed_data)
        # Adds warnings
        for message in parsed_data['warnings']:
            self.logger.error(message)
        
        # Parse dielectric function
        array_names_im = ['eps_energy', 'eps_im_x', 'eps_im_y', 'eps_im_z']  
        array_names_re = ['eps_energy', 'eps_re_x', 'eps_re_y', 'eps_re_z']  
        energy_units = 'eV'
        dielectric_units = 'arbitrary units'
        y_units = [dielectric_units]*3
        xy_data_im = XyData()
        xy_data_re = XyData()
                            
        try:
            dielectric_im_path = out_folder.get_abs_path(self._calc._OUTPUT_DIELECTRIC_IM)
            with open(dielectric_im_path, 'r') as fil:
                dielectric_im_file = fil.readlines()
            array_data = parse_raw_data(dielectric_im_file,array_names_im)
            xy_data_im.set_x(array_data[array_names_im[0]],array_names_im[0], energy_units) # energy
            y_arrays  = [array_data[array_names_im[1]],array_data[array_names_im[2]],array_data[array_names_im[3]]] # eps: x, y, z
            y_names = array_names_im[1:]
            xy_data_im.set_y(y_arrays,y_names,y_units)    
        except OSError:    
            self.logger.error("File with imaginary part of the dielectric function (interband) could not be found.")
            successful = False
            return successful, new_nodes_list
        try:
            dielectric_re_path = out_folder.get_abs_path(self._calc._OUTPUT_DIELECTRIC_RE)
            with open(dielectric_re_path, 'r') as fil:
                dielectric_re_file = fil.readlines()
            array_data = parse_raw_data(dielectric_re_file,array_names_re)
            xy_data_re.set_x(array_data[array_names_re[0]],array_names_re[0], energy_units) # energy
            y_arrays  = [array_data[array_names_re[1]],array_data[array_names_re[2]],array_data[array_names_re[3]]] # eps: x, y, z
            xy_data_re.set_y(y_arrays,y_names,y_units)    
        except OSError:    
            self.logger.error("File with real part of the dielectric function (interband) could not be found.")
            successful = False
            return successful, new_nodes_list
            
        new_nodes_list = [('output_parameters', output_params),
                          (self._eps_im_array_linkname, xy_data_im),
                          (self._eps_re_array_linkname , xy_data_re),
                          ]
        
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
        elif 'Fermi energy' in line and '=' in line:
            fermi_energy = line.split('=')[1].split('eV')[0]
            parsed_data['fermi_energy'] = fermi_energy
            parsed_data['fermi_energy_units'] = 'eV'
        elif 'Drude plasma frequency (xx)' in line:
            drude_plasma_freq_xx =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_xx'] = drude_plasma_freq_xx
            parsed_data['drude_plasma_frequency_units'] = 'eV'
        elif 'Drude plasma frequency (yy)' in line:
            drude_plasma_freq_yy =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_yy'] = drude_plasma_freq_yy
        elif 'Drude plasma frequency (zz)' in line:
            drude_plasma_freq_zz =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_zz'] = drude_plasma_freq_zz
        elif 'Drude plasma frequency (xy)' in line:
            drude_plasma_freq_xy =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_xy'] = drude_plasma_freq_xy
        elif 'Drude plasma frequency (xz)' in line:
            drude_plasma_freq_xz =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_xz'] = drude_plasma_freq_xz
        elif 'Drude plasma frequency (yz)' in line:
            drude_plasma_freq_yz =  line.split('=')[1].split('eV')[0]
            parsed_data['drude_plasma_frequency_yz'] = drude_plasma_freq_yz
            
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

def parse_raw_data(eps_file,array_names):
    """
    This function takes as input the eps_file as a list of filelines along
    with information on how to give labels and units to the parsed data
    
    :param dos_file: dos file lines in the form of a list
    :type dos_file: list
    :param array_names: list of all array names, note that array_names[0]
                        is for the case with non spin-polarized calculations
                        and array_names[1] is for the case with spin-polarized
                        calculation
    :type array_names: list
    :param array_units: list of all array units, note that array_units[0] is
                        for the case with non spin-polarized calculations and
                        array_units[1] is for the case with spin-polarized
                        calculation
    :type array_units: list
    
    :return array_data: narray, a dictionary for ArrayData type, which contains
                        all parsed dos output along with labels and units
    :return spin: boolean, indicates whether the parsed results are spin
                  polarized 
    """
    import numpy as np

    eps_header = eps_file[0]
    try:
        eps_data = np.genfromtxt(eps_file)
    except ValueError:
        raise QEOutputParsingError('epsfile could not be loaded '
        ' using genfromtxt')
    if len(eps_data) == 0:
        raise QEOutputParsingError("Dielectric file is empty.")
    if np.isnan(eps_data).any():
        raise QEOutputParsingError("Dielectric file contains non-numeric elements.")

    # Checks the number of columns, essentially to see whether spin was used
    if len(eps_data[0]) != 4:
        raise QEOutputParsingError("Eps file in format that the parser is not "
                                   "designed to handle.")

    i = 0
    array_data = {}
    array_data['header'] = np.array(eps_header)
    while i < len(array_names):
        array_data[array_names[i]] = eps_data[:, i]
        i += 1
    return array_data



