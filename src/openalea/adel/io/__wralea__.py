
# This file has been generated at Wed Mar 12 17:35:58 2014

from openalea.core import *


__name__ = 'openalea.adel.io'

__editable__ = True
__description__ = 'Manage input and output for adel such as meseasurement.'
__license__ = 'CECILL'
__url__ = ''
__alias__ = ['adel.io']
__version__ = '0.0.1'
__authors__ = 'C. Pradal, C. Fournier, C. Chambon'
__institutes__ = 'INRA, CIRAD, INRIA'
__icon__ = ''


__all__ = ['io_duplicate', 'io_canL2canS', 'io_mtg_factory', 'io_csvAsDict', 'io_RlistAsDict', 'io_dataframeAsdict', 'io_dataframe', 'io_thermal_time', 'io_to_canestra', 'io_saveRData', 'io_to_plantgl', 'PairAsDict_PairAsDict', 'io_readRData', 'io_apply_property', 'GetAdelString_GetAdelString', 'image_save_image', 'io_load_leaf_data', 'io_mtg2lpy', 'io_lpy2mtg', 'io_select_adel_geometric_data', 'io_select_adel_botanic_data']



io_duplicate = Factory(name='duplicate mtg',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='openalea.adel.io.io',
                nodeclass='duplicate',
                inputs=[{'name': 'mtg'}, {'interface': 'IInt', 'name': 'nb_plants', 'value': 1}],
                outputs=[{'name': 'mtg', 'desc': 'duplicated mtg'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_canL2canS = Factory(name='canL2canS',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='transform canopy Table Rdata to Surface canopy Table',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='canL2canS',
                inputs=[{'interface': None, 'name': 'canopy table (Rdataframe)', 'desc': ''}, {'interface': None, 'name': 'sr database (RList)', 'desc': ''}, {'interface': IFloat, 'name': 'senescence leaf shrink', 'value': 1, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_mtg_factory = Factory(name='mtg (params)',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='openalea.adel.io.io',
                nodeclass='mtg_factory',
                inputs=[{'interface': IDict, 'name': 'Canopy table'}, {'interface': IInt, 'name': 'number of sectors', 'value': 1}],
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




io_csvAsDict = Factory(name='CsvAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing the columns of the csv file',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='csvAsDict',
                inputs=[{'interface': IFileStr, 'name': 'Csv File name', 'desc': ''}, {'interface': IInt, 'name': 'csv type', 'value': 1, 'desc': '1 for coma separated/dot as decimal, 2 for ; separated, coma decimal'}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_RlistAsDict = Factory(name='RlistAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing the Robjects of the Rlist',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='RlistAsDict',
                inputs=[{'interface': IFileStr, 'name': 'R list', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_dataframeAsdict = Factory(name='dataframeAsdict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dictionary containing named vectors of values from adataframe (rpy2 object)',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='dataframeAsdict',
                inputs=[{'interface': IDict, 'name': 'dataframe', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_dataframe = Factory(name='Rdataframe',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dataframe (rpy2 object) from a dictionary containing named vectors of values',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='dataframe',
                inputs=[{'interface': IDict, 'name': 'Dict', 'desc': ''}],
                outputs=[{'interface': None, 'name': 'dataframe', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_thermal_time = Factory(name='update thermal time',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='simulation',
                nodemodule='openalea.adel.io.io',
                nodeclass='thermal_time',
                inputs=None,
                outputs=None,
                widgetmodule=None,
                widgetclass=None,
               )




io_to_canestra = Factory(name='to_canestra',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Convert Canestra Scene to a CanFile format stream.',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='to_canestra',
                inputs=[{'interface': IDict, 'name': 'scene', 'value': {}, 'desc': 'CanestraScene'}],
                outputs=[{'interface': ITextStr, 'name': 'string', 'desc': 'Canestra File'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_saveRData = Factory(name='saveRData',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='saveRData',
                inputs=[{'name': 'RObject', 'desc': 'R object'}, {'interface': IStr, 'name': 'Name of the savec object', 'value': 'Robj', 'desc': ''}, {'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'name', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_to_plantgl = Factory(name='to_plantgl',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Adapt a Canestra scene to a PlantGL scene',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='to_plantgl',
                inputs=[{'name': 'scene', 'desc': 'Canestra Scene'}, {'interface': IRGBColor, 'name': 'leaf_color', 'value': (0, 180, 0)}, {'interface': IRGBColor, 'name': 'stem_color', 'value': (0, 130, 0)}, {'interface': IRGBColor, 'name': 'soil_color', 'value': (170, 85, 0)}, {'interface': 'IDict', 'name': 'colors', 'desc': 'dict (vid, rgb color) '}, {'interface': 'IBool', 'name': 'ambient_only', 'value': False, 'desc': 'If True, set to 0 all optical properties except the ambient one'}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )




PairAsDict_PairAsDict = Factory(name='PairAsDict',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='returns a dict from a key,value pair',
                category='data i/o',
                nodemodule='openalea.adel.io.PairAsDict',
                nodeclass='PairAsDict',
                inputs=[{'interface': ITuple, 'name': 'Pair', 'value': None, 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_readRData = Factory(name='readRData',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='return a dictionary containing the Robject of the RData file',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='readRData',
                inputs=[{'interface': IFileStr, 'name': 'RData file', 'desc': ''}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_apply_property = Factory(name='apply_property',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='',
                category='data i/o, MTG',
                nodemodule='openalea.adel.io.io',
                nodeclass='apply_property',
                inputs=[{'interface': None, 'name': 'mtg', 'value': None, 'desc': 'MTG'}, {'interface': IStr, 'name': 'property name', 'value': '', 'desc': 'A MTG property'}, {'interface': IFunction, 'name': 'function', 'desc': 'function to apply on the property'}],
                outputs=[{'interface': IDict, 'name': 'dict', 'desc': 'the output property'}],
                widgetmodule=None,
                widgetclass=None,
               )




GetAdelString_GetAdelString = Factory(name='GetAdelString',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Reurn the content of a numbered lsystem string file',
                category='io',
                nodemodule='openalea.adel.io.GetAdelString',
                nodeclass='GetAdelString',
                inputs=[{'interface': IDirStr, 'name': 'LsysDir', 'value': None, 'desc': ''}, {'interface': IInt, 'name': 'IterNumber', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'LsysString', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




image_save_image = Factory(name='save_image',
                authors='C. Pradal, J. Coste ',
                description='Save a PlantGL scene in an image',
                category='io',
                nodemodule='openalea.adel.io.image',
                nodeclass='save_image',
                inputs=[{'name': 'scene', 'desc': 'PlantGL Scene'}, {'interface': 'IStr', 'name': 'image_name', 'value': '%s/img%04d.%s'}, {'interface': 'IDirStr', 'name': 'directory', 'value': '.'}, {'interface': 'IInt', 'name': 'index', 'value': 0}, {'interface': IEnumStr(enum=['png', 'jpg', 'tif']), 'name': 'ext', 'value': 'png'}],
                outputs=[{'interface': IInterface, 'name': 'scene', 'desc': 'PlantGL scene'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_load_leaf_data = Factory(name='load leaf data',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Load leaf data obtained by measurement',
                category='io',
                nodemodule='openalea.adel.io.io',
                nodeclass='load_leaf_data',
                inputs=[{'interface': IFileStr, 'name': 'filename', 'value': None, 'desc': 'filename of a pickle file'}],
                outputs=[{'interface': IDict, 'name': 'leaves', 'desc': 'A dictionary associating leaf rank to a list of leaf values (x, y, s, r)'}],
                widgetmodule=None,
                widgetclass=None,
               )




io_mtg2lpy = Factory(name='mtg2axial',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='Convert MTg to axial tree along with spec in lysystem',
                category='data i/o',
                nodemodule='openalea.adel.io.io',
                nodeclass='mtg2lpy',
                inputs=[{'interface': None, 'name': 'mtg', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'axial tree', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_lpy2mtg = Factory(name='lpy2mtg',
                authors='C. Pradal, C. Fournier (wralea authors)',
                description='aggregate lpy outputs into an mtg',
                category='data i/o',
                nodemodule='openalea.adel.io.io',
                nodeclass='lpy2mtg',
                inputs=[{'interface': None, 'name': 'axial tree', 'value': None, 'desc': ''}, {'interface': None, 'name': 'lsystem', 'value': None, 'desc': ''}, {'interface': None, 'name': 'scene', 'value': None, 'desc': ''}],
                outputs=[{'interface': None, 'name': 'mtg', 'desc': ''}],
                widgetmodule=None,
                widgetclass=None,
               )




io_select_adel_geometric_data = Factory(name='select_adel_geometric_data',
                authors='C.Chambon',
                description='',
                category='data i/o',
                nodemodule='openalea.adel.io.io',
                nodeclass='select_multiple_files',
                inputs=[{'interface': IStr, 'name': 'package', 'value': 'openalea.adel.data', 'desc': ''}, {'interface': IStr, 'name': 'laminaCur_pattern', 'value': '*laminaCur*.RData', 'desc': ''}, {'interface': IStr, 'name': 'lamina2D_pattern', 'value': '*lamina2D*.RData', 'desc': ''}, {'interface': IStr, 'name': 'laminaCur_filename', 'value': None, 'desc': ''}, {'interface': IStr, 'name': 'lamina2D_filename', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'laminaCur_filepath'}, {'interface': IStr, 'name': 'lamina2D_filepath'}],
                widgetmodule='io',
                widgetclass='SelectMultipleFiles',
               )




io_select_adel_botanic_data = Factory(name='select_adel_botanic_data',
                authors='C.Chambon',
                description='',
                category='data i/o',
                nodemodule='openalea.adel.io.io',
                nodeclass='select_multiple_files',
                inputs=[{'interface': IStr, 'name': 'package', 'value': 'openalea.adel.data', 'desc': ''}, 
                        {'interface': IStr, 'name': 'axeT_pattern', 'value': '*axeT*.csv', 'desc': ''}, 
                        {'interface': IStr, 'name': 'dimT_pattern', 'value': '*dimT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'phenT_pattern', 'value': '*phenT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'earT_pattern', 'value': '*earT*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'ssi2sen_pattern', 'value': '*ssi2sen*.csv', 'desc': ''},
                        {'interface': IStr, 'name': 'axeT_filename', 'value': None, 'desc': ''}, 
                        {'interface': IStr, 'name': 'dimT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'phenT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'earT_filename', 'value': None, 'desc': ''},
                        {'interface': IStr, 'name': 'ssi2sen_filename', 'value': None, 'desc': ''}],
                outputs=[{'interface': IStr, 'name': 'axeT_filepath'}, 
                         {'interface': IStr, 'name': 'dimT_filepath'},
                         {'interface': IStr, 'name': 'phenT_filepath'},
                         {'interface': IStr, 'name': 'earT_filepath'},
                         {'interface': IStr, 'name': 'ssi2sen_filepath'}],
                widgetmodule='io',
                widgetclass='SelectMultipleFiles',
               )

io_R_srdb = Factory(name='R_srdb',
                nodemodule='openalea.adel.AdelR',
                nodeclass='R_srdb',
               )
__all__.append('io_R_srdb')

io_R_xydb = Factory(name='R_xydb',
                nodemodule='openalea.adel.AdelR',
                nodeclass='R_xydb',
               )
__all__.append('io_R_xydb')
