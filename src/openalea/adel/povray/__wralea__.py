
# This file has been generated at Tue Oct  2 12:10:06 2012

from openalea.core import *


__name__ = 'openalea.adel.povray'

__editable__ = True
__description__ = 'Use of povray to compute coverage rate'
__license__ = 'CECILL'
__url__ = ''
__alias__ = []
__version__ = '0.9.1'
__authors__ = 'C. Pradal'
__institutes__ = 'CIRAD, INRIA'
__icon__ = ''


__all__ = ['povray_povray', 'povray_domain3D', 'povray_col_item', 'count_pixels_count_pixels', 'povray_stand_box']



povray_povray = Factory(name='povray',
                authors='C. Pradal (wralea authors)',
                description='',
                category='image',
                nodemodule='openalea.adel.povray.povray',
                nodeclass='povray',
                inputs=[{'name': 'scene', 'desc': 'PlantGL 3D scene'}, {'interface': IStr, 'name': 'pov_file', 'value': './scene.pov'}, {'interface': IFloat, 'name': 'camera_distance', 'value': 200}, {'interface': IInt, 'name': 'fov', 'value': 45}, {'interface': IInt, 'name': 'width', 'value': 320}, {'interface': IInt, 'name': 'height', 'value': 280}, {'interface': ITuple, 'name': 'domain', 'value': ((0, 0), (1, 1))}, {'interface': IInt, 'name': 'azimuth', 'value': 0}, {'interface': IInt, 'name': 'zenith', 'value': 0}, {'interface': IEnumStr(enum=['perspective', 'orthographic', 'fisheye']), 'name': 'camera_type', 'value': 'perspective'}, {'interface': IBool, 'name': 'soil', 'value': False}, {'interface': IStr, 'name': 'command', 'value': 'povray'}],
                outputs=[{'interface': IFileStr, 'name': 'povray image'}, {'interface': IFileStr, 'name': 'stand box image'}],
                widgetmodule=None,
                widgetclass=None,
               )




povray_domain3D = Factory(name='domain3D',
                authors='C. Pradal (wralea authors)',
                description='',
                category='image',
                nodemodule='openalea.adel.povray.povray',
                nodeclass='domain3D',
                inputs=[{'interface': 'ITuple', 'name': 'domain2D', 'value': ()}, {'interface': None, 'name': 'scene', 'value': None}],
                outputs=[{'interface': ITuple, 'name': 'domain3D'}],
                widgetmodule=None,
                widgetclass=None,
               )




povray_col_item = Factory(name='col_item',
                authors='C. Pradal (wralea authors)',
                description='',
                category='color',
                nodemodule='openalea.adel.povray.povray',
                nodeclass='col_item',
                inputs=[{'interface': 'IInt', 'name': 'color index', 'value': None, 'desc': 'color index. If None, return a function'}, {'interface': 'ISequence', 'name': 'color list', 'value': [(0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255), (255, 255, 0), (0, 255, 255), (255, 0, 255), (128, 255, 0), (0, 128, 255), (255, 0, 128), (0, 255, 128), (128, 0, 255), (255, 128, 0), (128, 128, 255), (255, 128, 128), (128, 255, 128), (255, 255, 255)]}],
                outputs=[{'interface': IRGBColor, 'name': 'color list'}],
                widgetmodule=None,
                widgetclass=None,
               )




count_pixels_count_pixels = Factory(name='count_pixels',
                authors='C. Chambon and B. Andrieu',
                description='',
                category='image processing',
                nodemodule='openalea.adel.povray.post_processing',
                nodeclass='count_pixels',
                widgetmodule=None,
                widgetclass=None,
               )




povray_stand_box = Factory(name='stand_box',
                authors='C. Pradal (wralea authors)',
                description='',
                category='image',
                nodemodule='povray',
                nodeclass='stand_box',
                inputs=[{'interface': 'ITuple', 'name': 'domain', 'value': ((0, 0, 0), (1, 1, 1))}],
                outputs=[{'interface': None, 'name': 'stand_box'}],
                widgetmodule=None,
                widgetclass=None,
               )




