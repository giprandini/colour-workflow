# -*- coding: utf-8 -*-
from aiida.orm.calculation.job.quantumespresso.namelists import NamelistsCalculation
from aiida.orm.calculation.job.quantumespresso.pw import PwCalculation
#TODO: remove imports, if they're really not being used
import os
from aiida.orm.calculation.job.quantumespresso import BasePwCpInputGenerator
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.common.utils import classproperty

__copyright__ = u"Copyright (c), This file is part of the AiiDA platform. For further information please visit http://www.aiida.net/. All rights reserved"
__license__ = "Non-Commercial, End-User Software License Agreement, see LICENSE.txt file."
__version__ = "0.7.0"
__authors__ = "The AiiDA team."

class SimpleCalculation(NamelistsCalculation):
    """
    Plugin for the simple.x code of the Quantum ESPRESSO distribution. 
    For more information regarding simple.x
    refer to ...
    """
    def _init_internal_params(self):
        super(SimpleCalculation, self)._init_internal_params()

        #self._SIMPLE_FILENAME = 'aiida.simple'
        self._default_namelists = ['INPUTSIMPLE']
        self._blocked_keywords = [#('SIMPLE','filsimple',self._SIMPLE_FILENAME),
                                  ('INPUTSIMPLE','outdir',self._OUTPUT_SUBFOLDER),
                                  ('INPUTSIMPLE','prefix',self._PREFIX),
                                 ]
        #self._internal_retrieve_list = [self._SIMPLE_FILENAME]
        self._default_parser = 'quantumespresso.simple'

    def use_parent_calculation(self,calc):
        """
        Set the parent calculation,
        from which it will inherit the outputsubfolder.
        The link will be created from parent RemoteData and NamelistCalculation
        """
        if not isinstance(calc,PwCalculation):
            raise ValueError("Parent calculation must be a PwCalculation")
        from aiida.common.exceptions import UniquenessError
        localdatas = [_[1] for _ in calc.get_outputs(also_labels=True)]
        if len(localdatas) == 0:
            raise UniquenessError("No output retrieved data found in the parent"
                                  "calc, probably it did not finish yet, "
                                  "or it crashed")

        localdata = [_[1] for _ in calc.get_outputs(also_labels=True)
                              if _[0] == 'remote_folder']
        localdata = localdata[0]
        self.use_parent_folder(localdata)
