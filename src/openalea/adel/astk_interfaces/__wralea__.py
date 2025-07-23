
# This file has been generated at Fri Nov 13 11:59:03 2015

from openalea.core import *


__name__ = 'openalea.adel.astk_interfaces'

__editable__ = True
__description__ = 'Interfaces complying to astk plant_interface'
__license__ = ''
__url__ = ''
__alias__ = []
__version__ = '0.0.1'
__authors__ = 'Christian Fournier'
__institutes__ = 'INRA'
__icon__ = ''


__all__ = ['openalea_adel_astk_interface_plot_statistics_node', 'openalea_adel_astk_interface_axis_statistics_node', 'openalea_adel_astk_interface_adelwheat_node', 'openalea_adel_geometric_elements_leaves_node', 'openalea_adel_astk_interface_adel_scene_node', 'openalea_adel_Stand_agronomicStand_node', 'openalea_adel_astk_interface_setup_canopy_node']



openalea_adel_astk_interface_plot_statistics_node = Factory(name='plot_statistics',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.astk_interface',
                nodeclass='plot_statistics_node',
                inputs=[{'interface': None, 'name': 'adel', 'value': None}, {'interface': None, 'name': 'Axis statistics', 'value': None}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_astk_interface_axis_statistics_node = Factory(name='axis_statistics',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.astk_interface',
                nodeclass='axis_statistics_node',
                inputs=[{'interface': None, 'name': 'adel', 'value': None}, {'interface': None, 'name': 'g', 'value': None}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_astk_interface_adelwheat_node = Factory(name='AdelWheat',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.astk_interface',
                nodeclass='adelwheat_node',
                inputs=[{'interface': IInt, 'name': 'nplants', 'value': 1}, {'interface': IInt, 'name': 'nsect', 'value': 1}, {'interface': None, 'name': 'devT', 'value': None}, {'interface': None, 'name': 'leaves', 'value': None}, {'interface': None, 'name': 'geoAxe', 'value': None}, {'interface': None, 'name': 'stand', 'value': None}, {'interface': None, 'name': 'run_adel_pars', 'value': None}, {'interface': IDict, 'name': 'options', 'value': {}}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_geometric_elements_leaves_node = Factory(name='Leaves',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.geometric_elements',
                nodeclass='leaves_node',
                inputs=[{'interface': None, 'name': 'xydb', 'value': None, 'desc': ''}, {'interface': None, 'name': 'srdb', 'value': None, 'desc': ''}, {'interface': None, 'name': 'geoLeaf', 'value': None, 'desc': ''}, {'interface': None, 'name': 'dynamic_bins', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'discretisation_level', 'value': 9, 'desc': ''}, {'interface': IInt, 'name': 'twist', 'value': 0, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'out', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_astk_interface_adel_scene_node = Factory(name='adel_scene',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.astk_interface',
                nodeclass='adel_scene_node',
                inputs=[{'interface': None, 'name': 'g', 'value': None}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_Stand_agronomicStand_node = Factory(name='AgronomicStand',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.Stand',
                nodeclass='agronomicStand_node',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




openalea_adel_astk_interface_setup_canopy_node = Factory(name='setup_canopy',
                authors='Christian Fournier (wralea authors)',
                description='',
                category='Unclassified',
                nodemodule='openalea.adel.astk_interface',
                nodeclass='setup_canopy_node',
                inputs=[{'interface': None, 'name': 'adel', 'value': None}, {'interface': IInt, 'name': 'age', 'value': 10}],
                outputs=({'interface': None, 'name': 'out'},),
                widgetmodule=None,
                widgetclass=None,
               )




