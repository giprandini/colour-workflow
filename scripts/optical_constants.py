# -*- coding: utf-8 -*-

from aiida.orm.calculation.inline import make_inline,optional_inline
from aiida.orm import DataFactory,Node,CalculationFactory
from aiida.backends.djsite.db import models
from aiida.orm.data.array import ArrayData
import numpy as np
import os , math
import json
from aiida.workflows.user.epfl_theos.dbimporters.utils import objects_set,objects_are_equal
from aiida.orm.utils import load_workflow
from tools import colour



__copyright__ = u"Copyright (c), This file is part of the AiiDA-EPFL Pro platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.1.0"
__authors__ = "Gianluca Prandini"


UpfData = DataFactory('upf')
ParameterData = DataFactory('parameter')


# Get dielectric function from ColourWorkflow
def prepare_input_optical_constants(workflow_pk):
    """
    workflow_pk: pk of ColourWorkflow
    """
    
    wf = load_workflow(workflow_pk)
    shirley_calc = wf.get_results()['shirley_calculation']
    drude_parameters = shirley_calc.out.output_parameters
    eps_im_inter = shirley_calc.out.array_eps_im
    eps_re_inter = shirley_calc.out.array_eps_re
    formula = wf.get_parameter('structure').get_formula()
    
    return drude_parameters, eps_im_inter, eps_re_inter, formula



# Write to file optical constants starting from the dielectric function calculated with the ColourWorkflow
def optical_constants_2file(parameters,drude_parameters,eps_im_inter,eps_re_inter,formula):
    """
    Calculate the optical constants
    :param parameters: ParameterData with input parameters
                        {'intra_broadening': [eV]}
    :param drude_parameters: output_parameters of a ShirleyCalculation 
    :param eps_im_inter: XyData from ShirleyCalculation  
    :param eps_re_inter: XyData from ShirleyCalculation
    :param formula: Name of the compound used in the files names
    :return: write to file optical constants
    """
    intra_broadening = parameters.get_dict()['intra_broadening']    
    
    drude_plasma_freq_x = float(drude_parameters.get_dict()['drude_plasma_frequency_xx'])
    drude_plasma_freq_y = float(drude_parameters.get_dict()['drude_plasma_frequency_yy'])
    drude_plasma_freq_z = float(drude_parameters.get_dict()['drude_plasma_frequency_zz'])
        
    energies = eps_im_inter.get_x()[1]  # eV  
    # eps interband (imaginary)
    eps_im_inter_x = eps_im_inter.get_y()[0][1]
    eps_im_inter_y = eps_im_inter.get_y()[1][1]
    eps_im_inter_z = eps_im_inter.get_y()[2][1]
    # eps interband (real)
    eps_re_inter_x = eps_re_inter.get_y()[0][1]
    eps_re_inter_y = eps_re_inter.get_y()[1][1]
    eps_re_inter_z = eps_re_inter.get_y()[2][1]
    # eps intraband
    eps_re_intra_x =  drude_plasma_freq_x**2/( energies**2 + intra_broadening**2 )
    eps_im_intra_x = drude_plasma_freq_x**2 * intra_broadening / ( energies * (energies**2 + intra_broadening**2)  )
    eps_re_intra_y =  drude_plasma_freq_y**2/( energies**2 + intra_broadening**2 )
    eps_im_intra_y = drude_plasma_freq_y**2 * intra_broadening / ( energies * (energies**2 + intra_broadening**2)  )
    eps_re_intra_z =  drude_plasma_freq_z**2/( energies**2 + intra_broadening**2 )
    eps_im_intra_z = drude_plasma_freq_z**2 * intra_broadening / ( energies * (energies**2 + intra_broadening**2)  )

    # eps interband+intraband
    eps_re_x = eps_re_inter_x - eps_re_intra_x
    eps_im_x = eps_im_inter_x + eps_im_intra_x
    eps_re_y = eps_re_inter_y - eps_re_intra_y
    eps_im_y = eps_im_inter_y + eps_im_intra_y
    eps_re_z = eps_re_inter_z - eps_re_intra_z
    eps_im_z = eps_im_inter_z + eps_im_intra_z
    # eps 
    eps_x = eps_re_x + 1j*eps_im_x
    eps_y = eps_re_y + 1j*eps_im_y
    eps_z = eps_re_z + 1j*eps_im_z
    
    #  EELS
    eels_x = eps_im_x / (eps_im_x**2 + eps_re_x**2)
    eels_y = eps_im_y / (eps_im_y**2 + eps_re_y**2)
    eels_z = eps_im_z / (eps_im_z**2 + eps_re_z**2)

    # eps average x,y,z
    eps = ( eps_x + eps_y + eps_z ) / 3.
    
    # The following quantities are all averaged along x,y,z
    norm_epsilon = np.sqrt(eps.real**2 + eps.imag**2)
    refractive_index = np.sqrt( ( eps.real + norm_epsilon  ) / 2. )
    extint_coeff = np.sqrt( ( -eps.real + norm_epsilon  ) / 2. )
    reflectivity = ( (refractive_index - 1.)**2 + extint_coeff**2 ) / ( (refractive_index + 1.)**2 + extint_coeff**2 )
    absorption = (energies*eps.imag/refractive_index)/1.9746*10.0e-7 #it is equivalent to 2.0*energies*extint_coeff/c or 4*pi*extint_coeff/lambda
    conductivity = 1j*energies/(4.*np.pi)*(1 - (eps.real + 1j*eps.imag) )   # conductivity = 1j*energies/(4.*np.pi)*(1 - eps)

    # Reflectivity as a function of the wavelength
    wavelengths_nm = 1239.8 / energies # in nm
    with open('reflectivity_{}_lambda.dat'.format(formula),'w') as o: 
        for i in xrange(len(wavelengths_nm)-1):
            o.write(str(wavelengths_nm[i]))
            o.write('  '+str(reflectivity[i]))
            o.write('\n')

    # Colours
    file_d65illuminant = '{}/D65_illuminant_1nm.dat'.format(os.path.dirname(colour.__file__))    
    file_cmf = '{}/cmf_1nm.dat'.format(os.path.dirname(colour.__file__)) 
    colours = colour.calcColour(energies, reflectivity, file_d65illuminant, file_cmf, do_plot=False)
    colours.pop('Fit_residuals')
   
    with open('colour_{}.dat'.format(formula),'w') as o:
        json.dump(colours, o, indent=4)            

    with open('drude_plasma_frequency_{}.dat'.format(formula),'w') as o:
        o.write('#     Drude plasma frequency (x,y,z) \n')
        o.write('{}   {}   {}'.format(str(drude_plasma_freq_x),str(drude_plasma_freq_y),str(drude_plasma_freq_z)))           

    with open('epsilon_inter_im_{}.dat'.format(formula),'w') as o: 
        o.write('#    Energy         x         y          z \n')
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(eps_im_inter_x[i]))
            o.write('  '+str(eps_im_inter_y[i]))
            o.write('  '+str(eps_im_inter_z[i]))
            o.write('\n')

    with open('epsilon_inter_re_{}.dat'.format(formula),'w') as o: 
        o.write('#    Energy         x         y          z \n')
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(eps_re_inter_x[i]))
            o.write('  '+str(eps_re_inter_y[i]))
            o.write('  '+str(eps_re_inter_z[i]))
            o.write('\n')

    with open('eels_{}.dat'.format(formula),'w') as o: 
        o.write('#    Energy         x         y          z \n')
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(eels_x[i]))
            o.write('  '+str(eels_y[i]))
            o.write('  '+str(eels_z[i]))
            o.write('\n')

    with open('epsilon_{}.dat'.format(formula),'w') as o: 
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(eps.imag[i]))
            o.write('  '+str(eps.real[i]))
            o.write('\n')

    with open('refractive_{}.dat'.format(formula),'w') as o: 
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(refractive_index[i]))
            o.write('  '+str(extint_coeff[i]))
            o.write('\n')

    with open('reflectivity_{}.dat'.format(formula),'w') as o: 
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(reflectivity[i]))
            o.write('\n')

    with open('conductivity_{}.dat'.format(formula),'w') as o: 
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(conductivity.imag[i]))
            o.write('  '+str(conductivity.real[i]))
            o.write('\n')

    # Drude plasma frequency averaged over the three Cartesian directions x, y, z
    drude_avg = math.sqrt( (drude_plasma_freq_x**2 + drude_plasma_freq_x**2 + drude_plasma_freq_z**2)/ 3. )
    with open('drude_plasma_frequency_avg_{}.dat'.format(formula),'w') as o:
        o.write('#     Drude plasma frequency (average) \n')
        o.write('{}'.format(str(drude_avg)))           

    # Dielectric function (interband contribution) averaged over the three Cartesian directions x, y, z
    eps_im_inter_avg = (eps_im_inter_x + eps_im_inter_y + eps_im_inter_z) / 3.
    eps_re_inter_avg = (eps_re_inter_x + eps_re_inter_y + eps_re_inter_z) / 3.
    with open('epsilon_inter_avg_{}.dat'.format(formula),'w') as o: 
        for i in xrange(len(energies)-1):
            o.write(str(energies[i]))
            o.write('  '+str(eps_im_inter_avg[i]))
            o.write('  '+str(eps_re_inter_avg[i]))
            o.write('\n')


if __name__ == "__main__":
    '''
    workflow_pk = pk of ColourWorkflow
    intra_broadening = empirical broadening for intraband contribution [eV]
    '''
    import sys , spglib
    workflow_pk = sys.argv[1]
    intra_broadening = sys.argv[2]
    (drude_parameters, eps_im_inter, eps_re_inter, formula) = prepare_input_optical_constants(int(workflow_pk))
    
    wf = load_workflow(int(workflow_pk))
    structure = wf.get_parameters()['structure']
    spacegroup_number = spglib.get_spacegroup(structure.get_ase()).split(' ')[1].split('(')[1].split(')')[0]

    os.mkdir('{}_pk{}_sg{}'.format(formula,workflow_pk,spacegroup_number))
    os.chdir('./{}_pk{}_sg{}'.format(formula,workflow_pk,spacegroup_number))

    parameters = ParameterData(dict={'intra_broadening': intra_broadening})  # specify value of gamma (intraband broadening in Drude-like expression)
    optical_constants_2file(parameters,drude_parameters,eps_im_inter,eps_re_inter,formula)


