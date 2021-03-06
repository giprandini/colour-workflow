#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from copy import deepcopy

'''
==================================================================================================
Input script to calculate the IPA dielectric function of elemental metals with the ColourWorkflow
==================================================================================================

sys.argv[1] = element

Notes:
- For aluminum we fix: interp_kmesh = 80x80x80, inter_broadening = 0.2 eV 

For more informations on the computational approach used and other technical details see: G. Prandini, PhD thesis, EPFL (2019)
'''

from time import sleep
import numpy as np
import sys , spglib
from aiida.common.example_helpers import test_and_get_code
from aiida.orm import  DataFactory
# from aiida.workflows.user.{username}.colour import ColourWorkflow
from aiida.workflows.user.prandini.colour import ColourWorkflow

StructureData = DataFactory('structure')
KpointsData = DataFactory('array.kpoints')
UpfData = DataFactory('upf')

def validate_upf_family(pseudo_family, all_species):
    """
    Validate if the pseudo_family is correct.
    """
    elements = all_species
    UpfData = DataFactory('upf')
    valid_families = UpfData.get_upf_groups(filter_elements=elements)
    valid_families_names = [family.name for family in valid_families]
    if pseudo_family not in valid_families_names:
        raise ValueError("Invalid pseudo family '{}'. "
                         "Valid family names are: {}".format(
            pseudo_family, ",".join(valid_families_names)))

def get_Zvalence_from_pseudo(pseudo):
     """
     Extract the number of valence electrons from a pseudopotential
     """
     with open(pseudo.get_file_abs_path(),'r') as f:
         lines=f.readlines()
         for line in lines:
             if 'valence' in line:
                 try:
                     return int(float(line.split("z_valence=\""
                                             )[-1].split("\"")[0].strip()))
                 except (ValueError, IndexError):
                     try:
                         return int(float(line.split("Z")[0].strip()))
                     except (ValueError, IndexError):
                         return None

## Dictionaries with wavefunction cutoffs for SG15 and PseudoDojo pseudopotentials
# Ecut for PseudoDojo (version 0.4). pseudo_family = 'Dojo_v0.4'
cutoffs_dojo_dict = {
'Au': 80.0,
'Rh': 95.0,
'K':  80.0,
'Ta': 65.0,
'Cu': 100.0,
'Ag': 90.0,
'Pd':90.0,
'Ti': 85.0,
'Ga': 85.0,
'Hf': 65.0,
'Mn': 100.0,
'Zn': 90.0,
'Fe': 100.0,
'Cr':100.0,
 }
# Ecut for SG15 (version 1.0 and 1.1). pseudo_family = 'SG15_1.0-1.1'
cutoffs_sg15_dict = {
'Al':80.0,
'Ag':55.0,
'As':50.0,
'Au':55.0,
'Be':60.0,
'Cd':55.0,
'Cu':90.0,
'Pd':55.0,
'Ga':120.0,
'W':55.0,
'K':55.0,
'Ca':55.0,
'Sc':55.0,
'Zr':55.0,
'Sn': 60.0,
'Y': 55.0,
'Sr': 55.0,
'Ti': 60.0,
'V': 70.0,
'Rh': 55.0,
'Na': 100.0,
'Li': 70.0,
'Nb': 90.0,
'Ta': 65.0,
'Re': 90.0,
'Mo': 50.0,
'Ir': 50.0,
'Pt': 65.0,
'Mg': 90.0,
'Rb': 50.0,
'In': 80.0,
'Sn': 65.0,
'Ru': 50.0,
'Tc': 55.0,
'Ni': 80.0,
'Fe': 80.0,  #  Delta factor is 4!
'Cr': 60.0,  #  Delta factor is 20!!!
'Si':50.0,
'Ge':70.0,
'Hf':70.0,
'Ba':60.0,
'C':80.0,
'N':80.0,
'P':65.0,
'Cs':60.0,
'B':55.0,  # Delta factor is 3 !
'Mn':65.0, # Delta factor is 13 !
'Sb':60.0,
'Te':55.0,
'Co':70.0,
'Os':70.0,
'Hg':65.0,
'Tl':60.0,
'Pb':50.0,
'Bi':50.0,
'Se':80.0,
}



############################################################################
## Input parameters
############################################################################

send = True

pw_codename = 'pw6.2@fidis'
simple_codename = 'simple6.2@fidis'
shirley_codename = 'simple_ip6.2@fidis'

element = sys.argv[1]

nscf_calculation_pk = 1259693
simple_calculation_pk = None

pseudo_family = 'SG15_1.0-1.1' 
#pseudo_family = 'Dojo_v0.4' 
dual = 4
relaxation_scheme = 'scf'
shirley_threshold = 0.0075
nonlocal_commutator = True   # for simple.x and simple_ip.x

# k-points
pw_kpoints_mesh = [24,24,24]   # scf mesh
nscf_kpoints_mesh = [2,2,2]    # nscf mesh
interp_kmesh = [64,64,64]      # interpolated mesh
if element == 'Al':
    interp_kmesh = [80,80,80]      # interpolated mesh


nscf_npools = 4           


###################################################################

# Crystal structures of the elemental metals are the Cottenier structures (K. Lejaeghere, Science 351, 2016) reduced to primitive cell
# AiiDA group with the crystal structures: 'Cottenier_structures_primitive'
if element in ['F']:
  raise Exception("Forbidden element")
g = Group.get_from_string('Cottenier_structures_primitive')
nodes = g.nodes.dbnodes
for node in nodes:
     name = node.attributes['kinds'][0]['name']
     if element == name:
        pk = node.pk
        break
     else:
        pass
structure = load_node(pk)

## validate
pw_code = test_and_get_code(pw_codename,'quantumespresso.pw')
test_and_get_code(simple_codename,'quantumespresso.simple')
test_and_get_code(shirley_codename,'shirley')
validate_upf_family(pseudo_family, structure.get_kind_names())

default_mpi_procs = pw_code.get_remote_computer().get_default_mpiprocs_per_machine()
# DFT resources
dft_max_num_machines = 2
dft_target_time_seconds = 60*60*1.0
dft_max_time_seconds = dft_target_time_seconds*1.0

# Simple resources
simple_max_num_machines = 2
simple_max_time_seconds = 60*60*4.0
simple_mpi_procs = default_mpi_procs 

# Shirley resources
shirley_max_num_machines = 8
shirley_max_time_seconds = 60*60*2.0
shirley_nthreads = 1   # Number of OpenMP threads
shirley_mpi_procs = 4  # default_mpi_procs / shirley_nthreads
omp_stacksize = '100m'

#if shirley_mpi_procs*shirley_nthreads != default_mpi_procs:
#    raise Exception("(MPI processes)*(OpenMP threads) is different from the number of CPUs")
if pseudo_family not in ['SG15_1.0-1.1','Dojo_v0.4']:
    raise Exception('Wrong pseudopotential family')

elements = structure.get_kind_names()
for element in elements:
  if element in ['Mn','Hf','Cr','Fe','Zn']:
    pseudo_family = 'Dojo_v0.4'
    validate_upf_family(pseudo_family, structure.get_kind_names())

spacegroup = spglib.get_spacegroup(structure.get_ase(), symprec=5e-4)
sg_number = spacegroup.split(' ')[1].split('(')[1].split(')')[0]
sg_name = spacegroup.split(' ')[0]
print 'Spacegroup={}'.format(spacegroup)

if nscf_calculation_pk:
   print 'Found previous nscf calculation (pk= {})'.format(nscf_calculation_pk)
if simple_calculation_pk:
   print 'Found previous simple calculation (pk= {})'.format(simple_calculation_pk)

# select wavefunction cutoff
if pseudo_family == 'SG15_1.0-1.1':
    wfc_cutoff = max([cutoffs_sg15_dict[element] for element in elements])
if pseudo_family == 'Dojo_v0.4':
    wfc_cutoff = max([cutoffs_dojo_dict[element] for element in elements])

print 'Wavefunction cutoff is {} Ry (pseudopotential family is {})'.format(wfc_cutoff,pseudo_family)
num_atoms_per_species = dict( (element,0) for element in elements )
for site in structure.get_site_kindnames():
   for element in elements:
      if site == element:
         num_atoms_per_species[element] += 1

# Find Z valence 
z_valence = {}
pseudo_family_group = UpfData.get_upf_group(pseudo_family)
for element in elements:
  for pseudo in pseudo_family_group.nodes:
    if element == pseudo.element:
       pseudo_element = pseudo
       break
    else:
       pass
  z_valence[element] = get_Zvalence_from_pseudo(pseudo_element)

number_of_electrons = 0.0
for element in elements:
   number_of_electrons += z_valence[element]*num_atoms_per_species[element]

# Number of bands (at least 30 empty bands)
occupied_bands = number_of_electrons/2.0
nbnd = int(occupied_bands) + max(30,int(number_of_electrons/1.333))
print 'nbnd={} (there are {} occupied bands)'.format(nbnd,int(occupied_bands))


# Workflow params
pw_calculation_set = {'custom_scheduler_commands': 'export OMP_NUM_THREADS=1'}
nscf_set_dict = {
              'custom_scheduler_commands': 'export OMP_NUM_THREADS=1',
              'resources':{'num_machines': dft_max_num_machines,
                            #'num_mpiprocs_per_machine':16
                            },
             "max_wallclock_seconds":dft_max_time_seconds,
             }
simple_set_dict = {
                   'custom_scheduler_commands': '#SBATCH --constraint=s6g1 \n#SBATCH  --mem=190000 \nexport OMP_NUM_THREADS=1',
                   'resources':{'num_machines': simple_max_num_machines,
                                'num_mpiprocs_per_machine':simple_mpi_procs,
                            },
                   'max_wallclock_seconds': simple_max_time_seconds,
                   }
shirley_set_dict = {'custom_scheduler_commands':
                     '#SBATCH --constraint=s6g1 \n#SBATCH --mem=190000 \nexport OMP_NUM_THREADS={}'.format(shirley_nthreads),
                    'resources':{'num_machines': shirley_max_num_machines,
                            'num_mpiprocs_per_machine':shirley_mpi_procs,
                            },
                    'max_wallclock_seconds':shirley_max_time_seconds,
                   }
shirley_settings = {'CMDLINE':['-c',str(shirley_nthreads)]}
#shirley_settings = {}

#npools = int(np.prod(nscf_kpoints_mesh)) 
#if npools > default_mpi_procs:
#    npools =  default_mpi_procs  # npools cannot be larger than the number of kpoints
nscf_settings = {'cmdline':['-nk',str(nscf_npools)]}
nscf_settings.update({ 'FORCE_KPOINTS_LIST':True, 'GAMMA_IMAGES':True})


# Pw input dictionary
pw_input_dict = {'CONTROL': {
                 'verbosity': 'high',
                 'tprnfor': True,
                 'tstress': True,
                 'disk_io': 'low',
                 },
              'SYSTEM': {
                 'ecutwfc': wfc_cutoff,
                 'ecutrho': wfc_cutoff*dual,
                 'occupations': 'smearing',
                 'degauss': 0.02,
                 'smearing': 'marzari-vanderbilt',
                 },
              'ELECTRONS': {
                 'mixing_beta': 0.3,
                 'scf_must_converge': True,   # not needed
                 },

                 }

nscf_input_dict = {'CONTROL': {
                             'verbosity': 'high',
                             'wf_collect': True,
                             },
                          'SYSTEM': {
                             'ecutwfc': wfc_cutoff,
                             'ecutrho': wfc_cutoff*dual,
                             'occupations': 'smearing',
                             'degauss': 0.02,
                             'smearing': 'marzari-vanderbilt',
                             'noinv': True,
                             'nosym': True,
                             'nbnd': nbnd,
                             },
                          'ELECTRONS': {
                             'mixing_beta': 0.3,
                             'diago_full_acc': True,
                             'diagonalization': 'cg',
                             },
                             }

simple_input_dict = {'INPUTSIMPLE':
                     {
                     's_bands':shirley_threshold,
                     'nkpoints(1)':interp_kmesh[0],
                     'nkpoints(2)':interp_kmesh[1],
                     'nkpoints(3)':interp_kmesh[2],
                     'nonlocal_commutator' : nonlocal_commutator,
                      }                     
                    }
shirley_input_dict = {'INPUTSIMPLEIP':
  {
   'simpleip_in%fermi_degauss' : 0.02205,
   'simpleip_in%fermi_ngauss' : -1,
   'simpleip_in%drude_degauss' : 0.00735,             # 0.1 eV
   'simpleip_in%wmin' : 0.0000735, 
   'simpleip_in%wmax' : 1.47,
   'simpleip_in%nw' : 2000,
   'simpleip_in%inter_broadening' : 0.00735,          # 0.1 eV	
   #'simpleip_in%inter_broadening' : 0.0147,          
   'simpleip_in%intra_broadening' : 0.00735,
   'simpleip_in%interp_grid(1)':interp_kmesh[0],
   'simpleip_in%interp_grid(2)':interp_kmesh[1],
   'simpleip_in%interp_grid(3)':interp_kmesh[2],
   'simpleip_in%nonlocal_commutator' : nonlocal_commutator,
  }
}

if element == 'Al':
    shirley_input_dict['INPUTSIMPLEIP']['simpleip_in%inter_broadening'] = 0.0147   # 0.2 eV (for Aluminum)

wf_parameters = {
                 'pw_codename': pw_codename,
                 'simple_codename': simple_codename,
                 'shirley_codename': shirley_codename,
                 'pseudo_family': pseudo_family,
                 'structure': structure,
                 'parameters':{
                               'pw_kpoints_mesh': pw_kpoints_mesh,
                               'nscf_kpoints_mesh': nscf_kpoints_mesh,
                               'dual': dual,
                               },
                 'input': { 
                           },

                 # Self-consistent Pw parameters
                 'pw_parameters': pw_input_dict,
                 'pw_calculation_set': pw_calculation_set, 
                 'pw_input':{'relaxation_scheme': relaxation_scheme,
                             'finish_with_scf': False,
#                              'volume_convergence_threshold': volume_conv_thr,
                             'automatic_parallelization':
                                                    {
                                                     'max_wall_time_seconds': dft_max_time_seconds,
                                                     'target_time_seconds': dft_target_time_seconds,
                                                     'max_num_machines': dft_max_num_machines
                                                     },
                          },                                
                 # nscf parameters
                 'nscf_parameters': nscf_input_dict,
#                  'nscf_interband_kpoints': nscf_interband_kpoints,
                'nscf_calculation_set': nscf_set_dict,
                'nscf_settings': nscf_settings,
                'nscf_input': {
                                          },
                 # simple parameters
                 'simple_parameters': simple_input_dict,
                 'simple_set_dict': simple_set_dict,
                 # shirley parameters
                 'shirley_parameters': shirley_input_dict,
                 'shirley_set_dict': shirley_set_dict,
                 'shirley_settings': shirley_settings,
                 }
if nscf_calculation_pk:
  wf_parameters.update({'nscf_calculation': load_node(nscf_calculation_pk)})
if simple_calculation_pk:
  wf_parameters.update({'simple_calculation': load_node(simple_calculation_pk)})
    

wf = ColourWorkflow(params = wf_parameters)

if send:
    wf.start()
    print ("Launch ColourWorkflow {}".format(wf.pk))
    print ("Parameters: {}".format(wf.get_parameters()))
else:
    print ("Would launch ColourWorkflow")
    print ("Parameters: {}".format(wf.get_parameters()))

if send:
    formula = structure.get_formula()
    with open('{}_{}_pk.shirley'.format(formula,relaxation_scheme),'a') as o:
        o.write(str(wf.pk))
        o.write('\n')
    sleep(1)   

