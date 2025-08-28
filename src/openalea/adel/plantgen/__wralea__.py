
# This file has been generated at Fri Oct 26 18:01:13 2012

from openalea.core import *


__name__ = 'openalea.adel.plantgen'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = ['plantgen']
__version__ = ''
__authors__ = 'M. Abichou, B. Andrieu, C. Chambon'
__institutes__ = ''
__icon__ = ''


__all__ = ['plantgen_plantgen', 'read_plantgen_inputs_read_plantgen_inputs', 'plantgen2adel_plantgen2adel', 'plot_HS_GL_SSI_T_plot_HS_GL_SSI_T', 'plot_dimT_plot_dimT', 'plot_tillering_dynamic_plot_tillering_dynamic', 'plot_tiller_probabilities_plot_tiller_probabilities']


plantgen_plantgen = Factory(name='plantgen',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='openalea.adel.plantgen.plantgen_interface',
                nodeclass='gen_adel_input_data',
                inputs=None,
                outputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'}, 
                         {'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'}, 
                         {'interface': None, 'name': 'phenT', 'desc': 'the phenT dataframe'}, 
                         {'interface': None, 'name': 'phenT_abs', 'desc': 'the phenT_abs dataframe'}, 
                         {'interface': None, 'name': 'dynT', 'desc': 'the dynT dataframe'}, 
                         {'interface': None, 'name': 'phenT_first', 'desc': 'the phenT_first dataframe'}, 
                         {'interface': None, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'}, 
                         {'interface': None, 'name': 'tilleringT', 'desc': 'the tilleringT dataframe'},
                         {'interface': None, 'name': 'cardinalityT', 'desc': 'the cardinalityT dataframe'},
                         {'interface': IDict, 'name': 'config', 'desc': 'the configuration used for the construction'}),
                widgetmodule=None,
                widgetclass=None,
               )

read_plantgen_inputs_read_plantgen_inputs = Factory(name='read_plantgen_inputs',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='openalea.adel.plantgen.plantgen_interface',
                nodeclass='read_plantgen_inputs',
                inputs=({'interface': IFileStr, 'name': 'inputs_filepath', 'desc': 'the Python module which describes the inputs of plantgen'},
                        {'interface': IFileStr, 'name': 'dynT_user_filepath', 'desc': 'the file path of the leaf dynamic parameters set by the user.'},
                        {'interface': IFileStr, 'name': 'dimT_user_filepath', 'desc': 'the file path of the dimensions of the organs set by the user.'}),
                outputs=({'interface': IInterface, 'name': 'dynT_user'},
                         {'interface': IInterface, 'name': 'dimT_user'},
                         {'interface': IInt, 'name': 'plants_number'},
                         {'interface': IInt, 'name': 'plants_density'},
                         {'interface': IDict, 'name': 'decide_child_axis_probabilities'},
                         {'interface': IDict, 'name': 'MS_leaves_number_probabilities'},
                         {'interface': IInt, 'name': 'ears_density'},
                         {'interface': IDict, 'name': 'GL_number'},
                         {'interface': IFloat, 'name': 'delais_TT_stop_del_axis'},
                         {'interface': IInterface, 'name': 'TT_hs_break'},
                         {'interface': IDict, 'name': 'inner_params'}),
                widgetmodule=None,
                widgetclass=None,
               )

plantgen2adel_plantgen2adel = Factory(name='plantgen2adel',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                category='data processing',
                nodemodule='openalea.adel.plantgen.plantgen_interface',
                nodeclass='plantgen2adel',
                inputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'},
                        {'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'},
                        {'interface': None, 'name': 'phenT', 'desc': 'the phenT dataframe'}),
                outputs=({'interface': None, 'name': 'axeT_adel', 'desc': 'the axeT dataframe in adel-like format'},
                         {'interface': None, 'name': 'dimT_adel', 'desc': 'the dimT dataframe in adel-like format'},
                         {'interface': None, 'name': 'phenT_adel', 'desc': 'the phenT dataframe in adel-like format'}),
                widgetmodule=None,
                widgetclass=None,
               )

plot_HS_GL_SSI_T_plot_HS_GL_SSI_T = Factory(name='plot_HS_GL_SSI_T',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                description='Plot the dimensions in `dimT` according to filters `dimensions_to_plot`, `id_dim_to_plot`, `id_cohort_to_plot`, `plot_non_regressive_tillers` and `plot_regressive_tillers`.',
                category='data processing',
                nodemodule='openalea.adel.plantgen.graphs',
                nodeclass='plot_HS_GL_SSI_T',
                inputs=({'interface': None, 'name': 'HS_GL_SSI_T', 'desc': 'the HS_GL_SSI_T dataframe'},
                        {'interface': None, 'name': 'dynT', 'desc': 'the dynT dataframe'},
                        {'interface': None, 'name': 'config', 'desc': 'a dictionary which stores the configuration used for the construction. `config` must contain at least the GL decimal numbers measured at several thermal times (including the senescence end).'},
                        {'interface': ISequence, 'name': 'id_phen_to_plot', 'value': [], 'desc': 'the list of id_phen to plot. If empty (the default) or None, then plot all the id_phen.'},
                        {'interface': ISequence, 'name': 'dynamics_to_plot',  'value': [], 'desc': "the list of dynamics to plot. The available dynamic are: 'HS', 'GL' and 'SSI'. If empty (the default) or None, then plot all the dynamics."},
                        {'interface': IBool, 'name': 'plot_most_frequent_main_stem', 'value': True, 'desc': 'whether to plot the most frequent main stem or not. If true, dynT must not be `None`.'},
                        {'interface': IBool, 'name': 'plot_non_regressive_tillers', 'value': True, 'desc': 'whether to plot the non regressive tillers or not. Non regressive tillers have id_dim ending by 1. Default is to plot the non regressive tillers.'},
                        {'interface': IBool, 'name': 'plot_regressive_tillers', 'value': True, 'desc': 'whether to plot the regressive tillers or not. Regressive tillers have id_dim ending by 0. Default is to plot the regressive tillers.'},
                        {'interface': IFileStr, 'name': 'plot_filepath', 'value': None, 'desc': 'the path of the file to save the plot in. If `None` (the default), do not save the plot but display it.'}),
                widgetmodule=None,
                widgetclass=None,
               )

plot_dimT_plot_dimT = Factory(name='plot_dimT',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                description='Plot the dimensions in `dimT` according to filters `dimensions_to_plot`, `id_dim_to_plot`, `id_cohort_to_plot`, `plot_non_regressive_tillers` and `plot_regressive_tillers`.',
                category='data processing',
                nodemodule='openalea.adel.plantgen.graphs',
                nodeclass='plot_dimT',
                inputs=({'interface': None, 'name': 'dimT', 'desc': 'the dimT dataframe'},
                        {'interface': ISequence, 'name': 'measured_id_dim', 'value': [], 'desc': 'the list of id_dim for which we have measured data. For these id_dim we plot no line but thick black markers. empty (the default) or None means that we do not represent measured data.'},
                        {'interface': IBool, 'name': 'relative_index_phytomer', 'value': False, 'desc': 'if True: display the relative index of the phytomers in abscissa. If False (the default), display the absolute index of the phytomers in abscissa.'},
                        {'interface': ISequence, 'name': 'dimensions_to_plot', 'value': [], 'desc': 'the list of dimensions to plot. If empty (the default) or None, then plot all the dimensions.'},
                        {'interface': ISequence, 'name': 'id_dim_to_plot', 'value': [], 'desc': 'the list of id_dim to plot. If empty (the default) or None, then plot all the id_dim.'},
                        {'interface': ISequence, 'name': 'id_cohort_to_plot', 'value': [], 'desc': 'the list of id_cohort to plot. If empty (the default) or None, then plot all the id_cohort.'},
                        {'interface': IBool, 'name': 'plot_non_regressive_tillers', 'value': True, 'desc': 'whether to plot the non regressive tillers or not. Non regressive tillers have id_dim ending by 1. Default is to plot the non regressive tillers.'},
                        {'interface': IBool, 'name': 'plot_regressive_tillers', 'value': True, 'desc': 'whether to plot the regressive tillers or not. Regressive tillers have id_dim ending by 0. Default is to plot the regressive tillers.'},
                        {'interface': IDirStr, 'name': 'plots_dirpath', 'value': None, 'desc': 'the path of the directory to save the plot(s) in.  If `None`, do not save the plot but display them.'}),
                widgetmodule=None,
                widgetclass=None,
               )


plot_tillering_dynamic_plot_tillering_dynamic = Factory(name='plot_tillering_dynamic',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                description='''Plot the dynamic of tillering, i.e. the evolution of the density of active axes when TT varies.
A non regressive axis is active at TT if TT >= TT_em_phytomer1. 
A regressive axis is active at TT if TT_em_phytomer1 <= TT < TT_stop_axis.''',
                category='data processing',
                nodemodule='openalea.adel.plantgen.graphs',
                nodeclass='plot_tillering_dynamic',
                inputs=({'interface': None, 'name': 'axeT', 'desc': 'the axeT dataframe'},
                        {'interface': IInt, 'name': 'plants_density', 'desc': 'the number of plants per square meter.'},
                        {'interface': IInt, 'name': 'TT_step', 'value': 10, 'desc': 'the thermal time step of the plot. Default is 10.'},
                        {'interface': IDirStr, 'name': 'plots_dirpath', 'value': None, 'desc': 'the path of the directory to save the plots in. If `None`, do not save the plots but display them.'}),
                widgetmodule=None,
                widgetclass=None,
               )


plot_tiller_probabilities_plot_tiller_probabilities = Factory(name='plot_tiller_probabilities',
                authors='M. Abichou, B. Andrieu, C. Chambon',
                description='''Plot the probability of each tiller.''',
                category='data processing',
                nodemodule='openalea.adel.plantgen.graphs',
                nodeclass='plot_tiller_probabilities',
                inputs=({'interface': None, 'name': 'dynT', 'desc': 'the dynT dataframe'},
                        {'interface': None, 'name': 'config', 'desc': 'a dictionary which stores the configuration used for the construction. `config` must contain at least the number of plants.'},
                        {'interface': IFileStr, 'name': 'plot_filepath', 'value': None, 'desc': 'the path of the file to save the plot in. If `None` (the default), do not save the plot but display it.'}),
                widgetmodule=None,
                widgetclass=None,
               )

