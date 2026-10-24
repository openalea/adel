
# This file has been generated at Wed Jul 29 18:55:10 2009

from openalea.core import *


__name__ = 'openalea.adel..lsystem'

__editable__ = True
__version__ = ''
__description__ = ''
__license__ = ''
__authors__ = ''
__url__ = ''
__institutes__ = ''


__all__ = ['_81622288', '_81621936', '_57070640', '_81622608', 'lsystem', '_81621648', '_82629808']


_81622288 = DataFactory(name='leaf.a',
                    description='',
                    editors={},
                    includes=None,
                    )


_81621936 = DataFactory(name='color.map',
                    description='',
                    editors={'Colormap': 'palette %s'},
                    includes=None,
                    )


_57070640 = DataFactory(name='description.txt',
                    description='',
                    editors={},
                    includes=None,
                    )


_81622608 = DataFactory(name='view.v',
                    description='',
                    editors={},
                    includes=None,
                    )



lsystem = CompositeNodeFactory(name='lsystem',
                             description='',
                             category='',
                             doc='',
                             inputs=[],
                             outputs=[],
                             elt_factory={  2: ('vlab.bin', 'process'),
   3: ('vlab.lsystem', 'Adel.l'),
   4: ('vlab.lsystem', 'description.txt'),
   5: ('vlab.lsystem', 'view.v'),
   6: ('vlab.lsystem', 'leaf.a'),
   7: ('vlab.lsystem', 'color.map'),
   8: ('vlab.lsystem', 'bid.txt')},
                             elt_connections={  9791872: (7, 0, 2, 0),
   9791884: (6, 0, 2, 0),
   9791896: (5, 0, 2, 0),
   9791908: (3, 0, 2, 0)},
                             elt_data={  2: {  'block': False,
         'caption': 'process',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([]),
         'posx': 250,
         'posy': 250,
         'priority': 0,
         'user_application': None},
   3: {  'block': False,
         'caption': 'Adel.l',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 125,
         'posy': 170,
         'priority': 0,
         'user_application': None},
   4: {  'block': False,
         'caption': 'description.txt',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 60,
         'posy': 40,
         'priority': 0,
         'user_application': None},
   5: {  'block': False,
         'caption': 'view.v',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 225,
         'posy': 170,
         'priority': 0,
         'user_application': None},
   6: {  'block': False,
         'caption': 'leaf.a',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 325,
         'posy': 170,
         'priority': 0,
         'user_application': None},
   7: {  'block': False,
         'caption': 'color.map',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 425,
         'posy': 170,
         'priority': 0,
         'user_application': None},
   8: {  'block': False,
         'caption': 'bid.txt',
         'hide': True,
         'lazy': True,
         'port_hide_changed': set([2]),
         'posx': 160,
         'posy': 40,
         'priority': 0,
         'user_application': None},
   '__in__': {  'block': False,
                'caption': 'In',
                'hide': True,
                'lazy': True,
                'port_hide_changed': set([]),
                'posx': 20,
                'posy': 5,
                'priority': 0,
                'user_application': None},
   '__out__': {  'block': False,
                 'caption': 'Out',
                 'hide': True,
                 'lazy': True,
                 'port_hide_changed': set([]),
                 'posx': 20,
                 'posy': 250,
                 'priority': 0,
                 'user_application': None}},
                             elt_value={  2: [(1, "'cpfg -m color.map Adel.l view.v leaf.a'")],
   3: [(0, 'PackageData(vlab.lsystem, Adel.l)'), (1, '{}'), (2, 'None')],
   4: [  (0, 'PackageData(vlab.lsystem, description.txt)'),
         (1, '{}'),
         (2, 'None')],
   5: [(0, 'PackageData(vlab.lsystem, view.v)'), (1, '{}'), (2, 'None')],
   6: [(0, 'PackageData(vlab.lsystem, leaf.a)'), (1, '{}'), (2, 'None')],
   7: [  (0, 'PackageData(vlab.lsystem, color.map)'),
         (1, "{'Colormap': 'palette %s'}"),
         (2, 'None')],
   8: [(0, 'PackageData(vlab.lsystem, bid.txt)'), (1, '{}'), (2, 'None')],
   '__in__': [],
   '__out__': []},
                             lazy=True,
                             )



_81621648 = DataFactory(name='bid.txt',
                    description='',
                    editors={},
                    includes=None,
                    )


_82629808 = DataFactory(name='Adel.l',
                    description='',
                    editors={},
                    includes=None,
                    )



