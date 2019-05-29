# -*- coding: utf-8 -*-
from __future__ import division
from copy import deepcopy
from aiida.common import aiidalogger
from aiida.orm.utils import load_node, load_workflow
from aiida.orm.workflow import Workflow
from aiida.orm import Code, WorkflowFactory, Group
from aiida.orm import CalculationFactory, DataFactory
from aiida.orm.data.array import ArrayData
from aiida.orm.implementation.django.calculation.job import JobCalculation
from aiida.workflows.user.epfl_theos.quantumespresso.pw import PwWorkflow, PwrestartWorkflow
from aiida.workflows.user.epfl_theos.quantumespresso.helpers import get_pw_wfs_with_parameters
from aiida.common.example_helpers import test_and_get_code

__copyright__ = u"Copyright (c), 2015, ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE (Theory and Simulation of Materials (THEOS) and National Centre for Computational Design and Discovery of Novel Materials (NCCR MARVEL)), Switzerland and ROBERT BOSCH LLC, USA. All rights reserved."
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.4.1"
__contributors__ = "Gianluca Prandini"

UpfData = DataFactory('upf')
ParameterData = DataFactory('parameter')
KpointsData = DataFactory('array.kpoints')
StructureData = DataFactory('structure')
PwCalculation = CalculationFactory('quantumespresso.pw')
SimpleCalculation = CalculationFactory('quantumespresso.simple')
logger = aiidalogger.getChild('ColourWorkflow')

## ===============================================
##    ColourWorkflow
## ===============================================

class ColourWorkflow(Workflow):
    """
    Workflow to compute the IPA dielectric function of a structure.
    
    Inputs:
    { 'pw_' + any parameter to be used in the pw workflow for the self-consistent DFT calculation
     
      'nscf_' + any parameter to be used in the pwrestart workflow for the nscf calculation
      
      'simple_' + any parameter to be used in the simple.x calculation

      'shirley_' + any parameter to be used in the simple_ip.x calculation
      
      'pseudo_family': pseudo family name,
      
      'parameters': {'pw_kpoints_mesh': list with k-point mesh for scf,
                     'nscf_kpoints_mesh': list with k-point mesh for nscf,
                     'dual': ecutrho / ecutwfc,
                    },
                
      'structure': AiiDA structure, 
      
     }
    
    """

    _default_volume_conv_thr = 0.01
    _default_conv_thr = 1.e-12
    
    def __init__(self, **kwargs):

        super(ColourWorkflow, self).__init__(**kwargs)


    ## ===============================================
    ##    Wf steps
    ## ===============================================


    @Workflow.step
    def start(self):
        """
        Check the input parameters
        """
        main_params = self.get_parameters()
        self.append_to_report("Starting the ColourWorkflow ...")
        self.append_to_report("Checking input parameters")

        try:
            main_params['nscf_calculation']
            has_nscf_calc = True
        except KeyError:
            has_nscf_calc = False
        try:
            main_params['simple_calculation']
            has_simple_calc = True
        except KeyError:
            has_simple_calc = False
        
        if has_nscf_calc and has_simple_calc:
            raise Exception('Found both nscf and simple calculation in the input of the workflow')

        
        self.next(self.run_pw)
        
    @Workflow.step
    def run_pw(self):
        """
        Run self-consistent DFT
        """
        main_params = self.get_parameters()
        pw_params = {}
        for k,v in main_params.iteritems():
            if k.startswith('pw_'):
                new_k = k[3:] # remove 'pw_' from the key name
                pw_params[new_k] = v
    
        pw_params['pseudo_family'] = main_params['pseudo_family']
        pw_params['codename'] = main_params['pw_codename']
        pw_params['structure'] = main_params['structure']
        if 'kpoints' not in pw_params:
            kpoints = KpointsData()
            kpoints.set_kpoints_mesh(main_params['parameters']['pw_kpoints_mesh'])
            kpoints.store()
            pw_params['kpoints'] = kpoints
        
        try:
            pw_params['input']['volume_convergence_threshold']
        except KeyError:
            pw_params['input']['volume_convergence_threshold'] = self._default_volume_conv_thr

        try:      
            pw_params['parameters']['ELECTRONS']['conv_thr']
        except KeyError:
            number_of_atoms = pw_params['structure'].get_ase().get_number_of_atoms()
            pw_params['parameters']['ELECTRONS']['conv_thr'] = self._default_conv_thr * number_of_atoms
        
        pw_wf_pks_already_computed = []  
        previous_pw_wfs = get_pw_wfs_with_parameters(pw_params, also_bands=False)
        previous_pw_wfs_states = [wf.get_state() for wf in previous_pw_wfs]
        if previous_pw_wfs_states and all(state == 'ERROR' for state in previous_pw_wfs_states):
                    self.append_to_report("The self-consistent calculation was already launched but it is in ERROR"
                                          " state. I'm launching again. Previous failed PwWorkflows pks"
                                          " are: {}".format([_.pk for _ in previous_pw_wfs]))
                    previous_pw_wfs = []
        if not previous_pw_wfs:
            wf = PwWorkflow(params=pw_params)
            self.attach_workflow(wf)
            wf.start()
            self.append_to_report("Launch PwWorkflow (pk: {}) for self-consistent"
                   " calculation".format(wf.pk))
            self.append_to_report("PwWorkflow Parameters: {}".format(wf.get_parameters()))
        else:
            if len(previous_pw_wfs)>1:
                self.append_to_report('WARNING! Found more than one previous PwWorkflow for the self-consistent'
                            ' calculation. I will take the last one.')
                previous_pw_wfs = [previous_pw_wfs[-1]]
            self.append_to_report("The self-consistent calculation was already done."
                   " Previous PwWorkflows pks are: {}".format([_.pk for _ in previous_pw_wfs]))
            pw_wf_pks_already_computed.extend([_.pk for _ in previous_pw_wfs if _.get_state() != 'ERROR'])
        self.add_attribute('pw_wf_pks_already_computed',pw_wf_pks_already_computed) 
        
        self.next(self.run_nscf)
    
    @Workflow.step    
    def run_nscf(self):    
        """
        Run non self-consistent DFT
        """  
        main_params = self.get_parameters()   

        try:
            main_params['nscf_calculation']
            has_nscf_calc = True
        except KeyError:
            has_nscf_calc = False
        try:
            main_params['simple_calculation']
            has_simple_calc = True
        except KeyError:
            has_simple_calc = False
        
        if not has_nscf_calc and not has_simple_calc:     
            list_of_wfs = [wf for wf in self.get_step(self.run_pw).get_sub_workflows()] \
                    + [load_workflow(pk) for pk in self.get_attribute('pw_wf_pks_already_computed')]
            if len(list_of_wfs)>1:
                raise Exception('Found more than one PwWorkflow for the previous self-consistent calculation') 
            if len(list_of_wfs)==0:
                raise Exception('Not found any previous PwWorkflow for the self-consistent calculation')
            pw_wf = list_of_wfs[0]
            pw_calc = pw_wf.get_result('pw_calculation')     
            
            nscf_wf_params = {}
            for k,v in main_params.iteritems():
                if k.startswith('nscf_'):
                    new_k = k[5:] # remove 'nscf_' from the key name
                    nscf_wf_params[new_k] = v
        
            nscf_wf_params['input']['relaxation_scheme'] = 'nscf'
            nscf_wf_params['pseudo_family'] = main_params['pseudo_family']
            nscf_wf_params['codename'] = main_params['pw_codename']
            nscf_wf_params['parameters']['CONTROL']['wf_collect'] = True     # write wfc to disk: necessary for simple.x
            nscf_wf_params['parameters']['ELECTRONS']['diago_full_acc'] = True 
            try:
                nscf_wf_params['structure'] = pw_wf.get_result('structure') 
            except ValueError:
                if main_params['pw_input']['relaxation_scheme'] == 'scf':
                    nscf_wf_params['structure'] = pw_calc.inp.structure 
                else:
                    raise Exception('Structure not found!')
            nscf_wf_params['remote_folder'] = pw_calc.out.remote_folder
            
            # k-points grid (uniform grid in the full 1BZ including the periodic images of Gamma)
            try:
                nscf_wf_params['settings'].update({'FORCE_KPOINTS_LIST':True, 'GAMMA_IMAGES':True})
            except KeyError:
                nscf_wf_params.update({'settings':{'FORCE_KPOINTS_LIST':True, 'GAMMA_IMAGES':True}})
            kpoints = KpointsData()
            kpoints.set_kpoints_mesh(main_params['parameters']['nscf_kpoints_mesh'])
            kpoints.store()
            nscf_wf_params['kpoints'] = kpoints
            
            nscf_wf_params['parameters']['SYSTEM']['nosym'] = True
            nscf_wf_params['parameters']['SYSTEM']['noinv'] = True
            
            try:      
                nscf_wf_params['parameters']['ELECTRONS']['conv_thr']
            except KeyError:
                number_of_atoms = pw_calc.res.number_of_atoms
                nscf_wf_params['parameters']['ELECTRONS']['conv_thr'] = self._default_conv_thr * number_of_atoms

            # If not specified in input_dict the default number of bands is put (arbitrarily we put at least 30 empty bands)    
            try:
                number_of_bands = nscf_wf_params['parameters']['SYSTEM']['nbnd']
            except KeyError:
                number_of_electrons = pw_calc.res.number_of_electrons
                number_of_bands = int(number_of_electrons/2) + max(30,int(0.3*number_of_electrons))
                nscf_wf_params['parameters']['SYSTEM']['nbnd'] = number_of_bands
        
            nscf_wf = PwrestartWorkflow(params=nscf_wf_params)
            self.attach_workflow(nscf_wf)
            nscf_wf.start()
            self.append_to_report("Launching PwWorkflow for 'nscf' calculation (pk: {})".format(nscf_wf.pk))
        elif has_nscf_calc and not has_simple_calc:
            self.append_to_report("Found previous nscf calculation (pk: {})".format(main_params['nscf_calculation'].pk))
        elif has_simple_calc:
            self.append_to_report("Found previous simple calculation: I will not do the nscf caclulation")    
        
                    
        self.next(self.run_simple)
            
    @Workflow.step    
    def run_simple(self):    
        """
        Run simple.x (construction of the optimal basis and computation of the relevant matrix elements)
        """
        main_params = self.get_parameters()

        try:
            main_params['simple_calculation']
            has_simple_calc = True
        except KeyError:
            has_simple_calc = False

        if has_simple_calc:
            self.append_to_report("Found previous simple calculation (pk: {})".format(main_params['simple_calculation'].pk))
        else:            
            try:
                nscf_calc = main_params.pop('nscf_calculation')
                if not isinstance(nscf_calc,PwCalculation):
                    raise Exception('nscf_calculation must be a QE calculation')
            except KeyError:
                nscf_wf = self.get_step_workflows(self.run_nscf)[0]
                nscf_calc = nscf_wf.get_result('pw_calculation')
                 
            if nscf_calc.inp.parameters.get_dict()['CONTROL']['calculation'] != 'nscf':
                raise Exception("The PwCalculation in wf_step 'run_nscf' must be a nscf calculation!")
                 
            num_of_electrons = int(nscf_calc.res.number_of_electrons)
            tot_num_of_bands = int(nscf_calc.res.number_of_bands)
            num_of_valence_bands = int(num_of_electrons / 2.)
            num_of_conduction_bands = tot_num_of_bands - num_of_valence_bands
            main_params['simple_parameters']['INPUTSIMPLE'].update({
                                                                    'calc_mode': 1,  # IP calculation
                                                                    'num_nbndv': num_of_valence_bands,
                                                                    'num_val': num_of_valence_bands,
                                                                    'num_cond': num_of_conduction_bands,
                                                                    })
            simple_code = Code.get_from_string(main_params['simple_codename'])
            simple_computer = simple_code.get_remote_computer()
            simple_parameters = ParameterData(dict = main_params['simple_parameters'])
            parentcalc = JobCalculation.get_subclass_from_pk(nscf_calc.pk)
            num_machines = main_params['simple_set_dict']['resources']['num_machines']
            max_wallclock_seconds = main_params['simple_set_dict']['max_wallclock_seconds']
            try:
                num_mpiprocs_per_machine = main_params['simple_set_dict']['resources']['num_mpiprocs_per_machine']
            except KeyError:
                pass
            try:
                custom_scheduler_commands = main_params['simple_set_dict']['custom_scheduler_commands']
            except KeyError:
                pass
             
            simple_calc = simple_code.new_calc(computer=simple_computer)
            simple_calc.set_max_wallclock_seconds(max_wallclock_seconds)
            try: 
                simple_calc.set_resources({"num_machines":num_machines,"num_mpiprocs_per_machine":num_mpiprocs_per_machine})
            except NameError:
                simple_calc.set_resources({"num_machines":num_machines})
            try:
                simple_calc.set_custom_scheduler_commands(custom_scheduler_commands)
            except NameError:
                pass

            simple_calc.use_parameters(simple_parameters)
            simple_calc.use_parent_calculation(parentcalc)
            simple_calc.store_all()
            self.append_to_report("Launching Simple calculation (pk: {})".format(simple_calc.pk))
            self.attach_calculation(simple_calc)

        self.next(self.run_shirley)
    
    @Workflow.step    
    def run_shirley(self):    
        """
        Run simple_ip.x for the calculation of the IPA dielectric function
        """
        main_params = self.get_parameters()

        try:
            simple_calc = main_params.pop('simple_calculation')
            if not isinstance(simple_calc,SimpleCalculation):
                raise Exception('simple_calculation is not a SimpleCalculation!')
        except KeyError:
            simple_calc = self.get_step_calculations(self.run_simple)[0]

        if simple_calc.get_state() == 'FINISHED':             
          shirley_code = Code.get_from_string(main_params['shirley_codename'])
          shirley_computer = shirley_code.get_remote_computer()
          shirley_parameters = ParameterData(dict = main_params['shirley_parameters'])
          parentcalc = JobCalculation.get_subclass_from_pk(simple_calc.pk)
          num_machines = main_params['shirley_set_dict']['resources']['num_machines']
          max_wallclock_seconds = main_params['shirley_set_dict']['max_wallclock_seconds']
          try:
            num_mpiprocs_per_machine = main_params['shirley_set_dict']['resources']['num_mpiprocs_per_machine']
          except KeyError:
            pass
          try:
            custom_scheduler_commands = main_params['shirley_set_dict']['custom_scheduler_commands']
          except KeyError:
            pass
          try:
            shirley_settings = main_params['shirley_settings']
          except KeyError:
            pass
         
          shirley_calc = shirley_code.new_calc(computer=shirley_computer)
          shirley_calc.set_max_wallclock_seconds(max_wallclock_seconds)
          try: 
            shirley_calc.set_resources({"num_machines":num_machines,"num_mpiprocs_per_machine":num_mpiprocs_per_machine})
          except NameError:
            shirley_calc.set_resources({"num_machines":num_machines})
          try:
            shirley_calc.set_custom_scheduler_commands(custom_scheduler_commands)
          except NameError:
            pass
          try:
            shirley_calc.use_settings(ParameterData(dict=shirley_settings))
          except NameError:
            pass
          shirley_calc.use_parameters(shirley_parameters)
          shirley_calc.use_parent_calculation(parentcalc)
          shirley_calc.store_all()
          self.append_to_report("Launching Shirley calculation (pk: {})".format(shirley_calc.pk))
          self.attach_calculation(shirley_calc)

          self.next(self.final_step)
        
        elif simple_calc.get_state() != 'FINISHED':
          self.append_to_report("WARNING! Shirley calculation (pk: {}) is in {} state. Stopping...".format(simple_calc.pk,simple_calc.get_state()))
          self.set_state('ERROR')
          self.next(self.exit)

    @Workflow.step    
    def final_step(self):    
        """
        
        """
        try:
            shirley_calc = self.get_step_calculations(self.run_shirley)[0]
        except IndexError:
            raise Exception('Shirley calculation was not performed!')
        # TODO: check that the calculation is finished correctly
        
        self.add_result('shirley_calculation', shirley_calc)
        
        drude_plasma_freq_x = float(shirley_calc.res.drude_plasma_frequency_xx)
        drude_plasma_freq_y = float(shirley_calc.res.drude_plasma_frequency_yy)
        drude_plasma_freq_z = float(shirley_calc.res.drude_plasma_frequency_zz)
        
        array_eps_im = shirley_calc.out.array_eps_im
        array_eps_re = shirley_calc.out.array_eps_re
        if set(array_eps_im.get_x()[1] == array_eps_re.get_x()[1]) != {True}:
            raise Exception('Inconsistency in the energy grid of the dielectric function')
        
        array_eps_im.label = "Imaginary part of the interband dielectric function along x,y and z"
        array_eps_im.description = ("Calculated with"
                                 " the ColourWorkflow (pk: {})".format(self.pk))
        array_eps_re.label = "Real part of the interband dielectric function along x, y and z"
        array_eps_re.description = ("Calculated with"
                                 " the ColourWorkflow (pk: {})".format(self.pk))
        self.append_to_report("Saved imaginary part of the interband dielectric function (XyData pk: {})"
                                  "".format(array_eps_im.pk))
        self.add_result("array_eps_im", array_eps_im)
        self.append_to_report("Saved real part of the interband dielectric function (XyData pk: {})"
                                  "".format(array_eps_re.pk))
        self.add_result("array_eps_re", array_eps_re)
        
        self.add_result('drude_plasma_frequency_x', drude_plasma_freq_x)
        self.add_result('drude_plasma_frequency_y', drude_plasma_freq_y)
        self.add_result('drude_plasma_frequency_z', drude_plasma_freq_z)
        
        self.append_to_report("Ending the ColourWorkflow.")
        self.next(self.exit)



   
