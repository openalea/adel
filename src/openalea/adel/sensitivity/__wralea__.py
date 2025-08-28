
# This file has been generated at Thu Apr 30 14:12:27 2009

from openalea.core import *


__name__ = 'openalea.adel.sensitivity'

__editable__ = True
__description__ = 'Utilities for sensitivity analysis'
__license__ = 'CECILL'
__url__ = ''
__version__ = '0.0.1'
__authors__ = 'C. Fournier, C. Pradal'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['sensitivity_Morris', 'sensitivity_MorrisIS','sensitivity_plotSens']



sensitivity_Morris = Factory(name='Morris',
                description='return a simulation plan for Morris sensitivity analysis',
                category='',
                nodemodule='openalea.adel.sensitivity.sensitivity',
                nodeclass='Morris',
                inputs=[{'interface': IInt, 'name': 'NbTrajectories', 'value': 10, 'desc': 'Number of trajectory (nb run = (nb factors+1)*nbTrajectories'},
                        {'interface': ISequence, 'name': 'factors', 'desc': 'list of naes of factors to vary'},
                        {'interface': ISequence, 'name': 'binf', 'desc': 'inferior bound of factors'},
                        {'interface': ISequence, 'name': 'bsup', 'desc': 'Superior bound of factors'}],
                outputs=[{'interface': None, 'name': 'MorrisPlan', 'desc': 'Simulated plan (Robject)'},
                         {'interface': None, 'name': 'Parameters', 'desc': 'Parameters'}],
                widgetmodule=None,
                widgetclass=None,
               )


sensitivity_MorrisIS = Factory(name='MorrisIS',
                description='Compute sensitivity indices for a Morris simulation plan',
                category='',
                nodemodule='openalea.adel.sensitivity.sensitivity',
                nodeclass='Morris_IS',
                inputs=[{'interface': None, 'name': 'MorrisPlan', 'desc': 'Simulated plan (Robject)'},
                        {'interface': ISequence, 'name': 'model response', 'desc': 'list of Y (model output)'}],
                outputs=[{'interface': None, 'name': 'MorrisOut', 'desc': 'Result of sensitivity (Robject)'},
                         {'interface': None, 'name': 'mustar', 'desc': 'mean of absolute factors effect'},
                         {'interface': None, 'name': 'sigma', 'desc': 'standard error of factors effect'}],
                widgetmodule=None,
                widgetclass=None,
               )

sensitivity_plotSens = Factory(name='plotSens',
                description='plot sensitivity analysis object (R)',
                category='',
                nodemodule='openalea.adel.sensitivity.sensitivity',
                nodeclass='plotSens',
                inputs=[{'interface': None, 'name': 'SensOut', 'desc': 'Simulated sensitivity Robject'}],
                outputs=[{'interface': None, 'name': 'out'}],
                widgetmodule=None,
                widgetclass=None,
               )


