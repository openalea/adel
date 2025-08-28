
# This file has been generated at Tue Feb 24 14:26:36 2015

from openalea.core import *


__name__ = 'openalea.adel.tutorials.ssdbm.leaves'

__editable__ = True
__description__ = ''
__license__ = ''
__url__ = ''
__alias__ = ['ssdbm.leaves']
__version__ = ''
__authors__ = 'C. Pradal, C. Fournier, S. Cohen-Boulakia'
__institutes__ = 'INRIA'
__icon__ = ''


__all__ = ['leaves']



leaves = CompositeNodeFactory(name='leaves',
                             description='Fit leaf data into a geometrical model',
                             category='geometry',
                             doc='',
                             inputs=[  
                             {'desc': 'Midrib curves of leaves', 'interface': IFileStr, 'name': 'leaf_database', 'value': None},
                             {'desc': 'Number of leaf sectors', 'interface': IInt, 'name': 'nb_sectors', 'value': 7}],
                             outputs=[{  'desc': '', 'interface': None, 'name': 'symbols'}],
                             elt_factory={  2: ('openalea.adel.io', 'load leaf data'),
   3: ('openalea.adel.fitting', 'fit leaves'),
   4: ('openalea.adel.geometry', 'symbols')},
                             elt_connections={  140481681999728: ('__in__', 1, 3, 1),
   140481681999752: ('__in__', 0, 2, 0),
   140481681999776: (4, 0, '__out__', 0),
   140481681999800: (2, 0, 3, 0),
   140481681999824: (3, 0, 4, 0)},
                             elt_data={  2: {  'block': False,
         'caption': 'load leaf data',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x37c8c50> : "load leaf data"',
         'hide': True,
         'id': 2,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': -166.74502103784138,
         'posy': -109.95444946135717,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   3: {  'block': False,
         'caption': 'fit leaves',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x5488650> : "fit leaves"',
         'hide': True,
         'id': 3,
         'is_in_error_state': False,
         'is_user_application': False,
         'lazy': True,
         'minimal': False,
         'port_hide_changed': set([]),
         'posx': -92.4227816421706,
         'posy': -43.20169118234779,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   4: {  'block': False,
         'caption': 'symbols',
         'delay': 0,
         'factory': '<openalea.core.node.NodeFactory object at 0x6991b10> : "symbols"',
         'hide': True,
         'id': 4,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': -109.65382073133082,
         'posy': 41.21903888954781,
         'priority': 0,
         'use_user_color': False,
         'user_application': None,
         'user_color': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'delay': 0,
                'hide': True,
                'id': 0,
                'lazy': True,
                'port_hide_changed': set([]),
                'posx': -74.9879298828165,
                'posy': -169.97674397328998,
                'priority': 0,
                'use_user_color': False,
                'user_application': None,
                'user_color': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'delay': 0,
                 'hide': True,
                 'id': 1,
                 'lazy': True,
                 'port_hide_changed': set([]),
                 'posx': -96.45061824038132,
                 'posy': 72.99512198066596,
                 'priority': 0,
                 'use_user_color': False,
                 'user_application': None,
                 'user_color': None}},
                             elt_value={  2: [], 3: [], 4: [(1, 'None'), (2, 'False')], '__in__': [], '__out__': []},
                             elt_ad_hoc={  2: {'useUserColor': False, 'position': [-166.74502103784138, -109.95444946135717], 'userColor': None},
   3: {'useUserColor': False, 'position': [-92.4227816421706, -43.20169118234779], 'userColor': None},
   4: {'useUserColor': False, 'position': [-109.65382073133082, 41.21903888954781], 'userColor': None},
   '__in__': {'useUserColor': False, 'position': [-74.9879298828165, -169.97674397328998], 'userColor': None},
   '__out__': {'useUserColor': False, 'position': [-96.45061824038132, 72.99512198066596], 'userColor': None}},
                             lazy=True,
                             eval_algo='LambdaEvaluation',
                             )




