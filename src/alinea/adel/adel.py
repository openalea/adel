"""Mother class for Adel models"""
import os
import pandas
import numpy
import warnings
import operator
from itertools import chain

try:
    import cPickle as pickle
except:
    import pickle
from openalea.plantgl.all import Viewer, Scene
from openalea.mtg.plantframe.color import colormap

from alinea.adel.geometric_elements import Leaves
from alinea.adel.Stand import AgronomicStand
from alinea.adel.mtg_interpreter import plot3d, transform_geom, mtg_interpreter
from alinea.adel.postprocessing import axis_statistics, plot_statistics, \
    midrib_statistics
from alinea.adel.newmtg import exposed_areas, exposed_areas2canS, duplicate, \
    mtg_factory


def flat_list(nested_list):
    return list(chain.from_iterable(nested_list))


def balanced_sample(n, proba):
    """ return a list of n keys found in proba, repecting probalities of proba values"""
    card = {k: int(v * n) for k, v in proba.items()}
    missing = int(n - sum(card.values()))
    while (missing > 0):
        # current frequencies
        freq = {k: float(v) / n for k, v in card.items()}
        # diff with probabilities
        dp = {k: abs(freq[k] - proba[k]) for k in freq}
        sorted_p = sorted(dp.items(), key=operator.itemgetter(1), reverse=True)
        k = sorted_p[0][0]
        card[k] += 1
        missing -= 1
    card = {k: v for k, v in card.items() if v > 0}
    items = flat_list([[key] * val for key, val in card.items()])
    numpy.random.shuffle(items)
    return items


class Adel(object):
    """Mother class for adel models"""
    conv_units = {'mm': 0.001, 'cm': 0.01, 'dm': 0.1, 'm': 1, 'dam': 10,
                  'hm': 100,
                  'km': 1000}

    def __init__(self, nref_plants=1, nplants=1, duplicate=None, species=None,
                 nsect=1,
                 leaves=None, stand=None,
                 aspect='smart', split=False,
                 face_up=False, classic=False, scene_unit='cm',
                 age=None, seed=None, leaf_db=None, positions=None,
                 convUnit=None):
        """

        Args:
            nref_plants: the number of reference plants in the database
            nplants: the number of plants in the canopy
            duplicate:
            species: a {species: frequency} dict indicating the composition of
             the canopy. If None (default), a monospecific canopy of species '0'
             is generated
            nsect: (int) the number of sectors on leaves
            leaves: (object) a Leaves class instance pointing to leaf shape
             database or a {species:leaf_db} dict referencing distinct database
             per species
            stand: (object) a Stand class instance
            aspect: (str) the aspect of the stand (square, line or smart)
            split:
            face_up:
            classic: (bool) should stem cylinders be classical pgl cylinders ?
            scene_unit: (string) desired length unit for the output mtg
            age: (optional) the age of the canopy
            seed: (int) a seed for the random number generator
            leaf_db: deprecated, use leaves
            positions: deprecated, use stand
            convUnit: deprecated, use scene_unit
        """

        if leaf_db is not None:
            warnings.warn(
                '!!!!Warning!!!! leaf_db argument is deprecated, '
                'use adel.geometric_elements.Leaves class instead',
                DeprecationWarning)

        if positions is not None:
            warnings.warn(
                '!!!!Warning!!!! positions argument is deprecated,'
                ' use stand = adel.Stand class instead',
                DeprecationWarning)

        if convUnit is not None:
            warnings.warn(
                "!!!!Warning!!!! convUnit argument is deprecated,"
                " use scene_unit=unit_name (" + ','.join(self.conv_units) +
                ") instead",
                DeprecationWarning)

        if species is None:
            species = {0: 1}

        if leaves is None:
            leaves = {0: Leaves()}

        if not isinstance(leaves, dict):
            leaves = {0: leaves}

        if stand is None:
            stand = AgronomicStand(sowing_density=250, plant_density=250,
                                   inter_row=0.15)

        self.nref_plants = nref_plants
        self.nplants = nplants
        self.duplicate = duplicate
        self.stand = stand
        self.species = species
        self.leaves = leaves
        self.scene_unit = scene_unit
        self.convUnit = self.conv_units[self.scene_unit]
        self.nsect = nsect
        self.split = split
        self.face_up = face_up
        self.classic = classic
        self.seed = seed
        self.meta = {}

        self.new_stand(nref_plants=nref_plants, nplants=nplants,
                       duplicate=duplicate, seed=seed,
                       aspect=aspect, age=age, species=species)

    def new_stand(self, nref_plants=None, nplants=None, duplicate=None,
                  seed=None, aspect=None,
                  age=None, species=None):

        if nref_plants is not None:
            self.nref_plants = nref_plants

        if nplants is not None:
            self.nplants = nplants

        if seed is not None:
            self.seed = seed
            numpy.random.seed(self.seed)

        if age is None:
            self.canopy_age = -999
        else:
            self.canopy_age = age
        self.meta.update({'canopy_age': self.canopy_age})

        if species is not None:
            self.species = species

        if aspect is not None:
            self.aspect = aspect

        if duplicate is not None:
            self.duplicate = duplicate

        if self.aspect is 'smart':
            self.nplants, self.domain, self.positions, \
            self.domain_area = self.stand.smart_stand(
                self.nplants, at=age, convunit=1. / self.convUnit)
        else:
            self.nplants, self.domain, self.positions, \
            self.domain_area = self.stand.stand(
                self.nplants, aspect=self.aspect, convunit=1. / self.convUnit)

        self.plant_azimuths = numpy.random.random(self.nplants) * 360
        self.plant_species = balanced_sample(self.nplants, self.species)
        self.plant_references = numpy.random.choice(range(self.nref_plants),
                                                    self.nplants)

        stand_parameters = {'sowing_density': self.stand.sowing_density,
                            'plant_density': self.stand.plant_density,
                            'inter_row': self.stand.inter_row,
                            'noise': self.stand.noise,
                            'density_curve_data': self.stand.density_curve_data}
        self.meta.update(
            {'stand': stand_parameters, 'aspect': self.aspect,
             'nplants': self.nplants,
             'domain': self.domain, 'domain_area': self.domain_area,
             'nsect': self.nsect, 'scene_unit': self.scene_unit,
             'convUnit': self.convUnit,
             'split': self.split})

        if duplicate is not None:
            # split nplants into nquot and nrem, so that
            # nplants = nquote * duplicate + nrem
            self.nrem = self.nplants % duplicate
            self.nquot = (self.nplants - self.nrem) / duplicate
            if self.nquot == 0 and self.duplicate > 0:  # degenerated case
                self.nquot = 1
                self.nrem = 0
                self.duplicate = self.nplants

    def duplicated(self, gquot, grem=None):
        """Construct g using duplications"""
        if self.duplicate is None:
            raise ValueError('Duplication not defined for this stand')
        g = duplicate(gquot, self.nquot * self.duplicate, grem)
        # dispose plants and renumber them
        pos = g.property('position ')
        az = g.property('azimuth')
        lab = g.property('label')
        geom = g.property('geometry')
        for i, vid in enumerate(g.vertices(1)):
            lab[vid] = 'plant' + str(i + 1)
            pos[vid] = self.positions[i]
            az[vid] = self.plant_azimuths[i]
            for gid in g.components_at_scale(vid, g.max_scale()):
                if gid in geom:
                    geom[gid] = transform_geom(geom[gid], self.positions[i],
                                               self.plant_azimuths[i])
        return g


    def build_mtg(self, parameters, stand, **kwds):
        g = mtg_factory(parameters, stand=stand, leaf_sectors=self.nsect,
                        leaves=self.leaves, split=self.split, **kwds)
        g = mtg_interpreter(g, self.leaves, classic=self.classic,
                            face_up=self.face_up)
        return g


    def meta_informations(self, g):
        if 'meta' in g.property_names():
            return g.property('meta').values()[0]
        else:
            return self.meta

    @staticmethod
    def get_axis(g, plant='plant1', axe='MS'):
        """ return a new mtg representing an axe
        """
        p = [vid for vid in g.vertices(scale=1) if g.label(vid) == plant][0]
        ax = [vid for vid in g.components(p) if g.label(vid) == axe][0]
        return g.sub_mtg(ax, copy=True)

    @staticmethod
    def plot(g, property=None):
        s = Adel.scene(g, property)
        Viewer.display(s)
        return s

    @staticmethod
    def scene(g, property=None):
        if property:
            g = colormap(g, property, cmap='jet', lognorm=True)
            colored = g.property('color')
            colors = {
                vid: colored.get(vid, colored.get(g.complex(vid), [0, 0, 0]))
                for
                vid in g.vertices(scale=g.max_scale())}
        else:
            colors = None
        s = plot3d(g,
                   colors=colors)  # use the one of openalea.plantframe.color instead ?
        return s

    def get_exposed_areas(self, g, convert=False, TT=None):
        areas = exposed_areas(g)
        if convert:
            areas = exposed_areas2canS(areas)
        if TT is None:
            TT = self.meta_informations(g)['canopy_age']
        areas['TT'] = TT
        return areas

    def axis_statistics(self, g):
        meta = self.meta_informations(g)
        df_lai = self.get_exposed_areas(g, convert=True)
        axstat = None
        if not df_lai.empty:
            axstat, _ = axis_statistics(df_lai, meta['domain_area'],
                                        meta['convUnit'])
        return axstat

    def plot_statistics(self, g=None, axstat=None):
        meta = self.meta_informations(g)
        pstat = None
        if axstat is None:
            axstat = self.axis_statistics(g)
        pstat = plot_statistics(axstat, meta['nplants'], meta['domain_area'])
        return pstat

    def save(self, g, index=0, dir='./adel_saved', basename=None,
             check_meta=True):
        if check_meta:
            if 'meta' not in g.property_names():
                root = g.node(0)
                root.meta = self.meta
        if basename is None:
            if not os.path.exists(dir):
                os.mkdir(dir)
            basename_geom = dir + '/scene%04d' % (index)
            basename_adel = dir + '/adel%04d' % (index)
        else:
            basename_adel = basename_geom = str(basename)
        s = Adel.scene(g)
        geom = {sh.id: sh.geometry for sh in s}
        g.remove_property('geometry')
        fgeom = basename_geom + '.bgeom'
        fg = basename_adel + '.pckl'
        s.save(fgeom, 'BGEOM')
        with open(fg, 'wb') as output:
            pickle.dump(g, output, protocol=2)
        # restore geometry
        g.add_property('geometry')
        g.property('geometry').update(geom)
        return fgeom, fg

    @staticmethod
    def load(index=0, directory='./adel_saved', basename=None, load_geom=True):
        if basename is None:
            if not os.path.exists(dir):
                os.mkdir(dir)
            basename_geom = dir + '/scene%04d' % (index)
            basename_adel = dir + '/adel%04d' % (index)
        else:
            basename_adel = basename_geom = basename
        fgeom = basename_geom + '.bgeom'
        fg = basename_adel + '.pckl'
        if not os.path.exists(fgeom) or not os.path.exists(fg):
            raise IOError('adel cannot find saved files')

        f = open(fg, 'rb')
        g = pickle.load(f, encoding='bytes')
        f.close()

        # backward compatibility
        if isinstance(g, list):
            g, age = g
            root = g.node(0)
            meta = {'canopy_age': age}
            if 'meta' not in g.property_names():
                root.meta = meta
            else:
                root.meta.update(meta)

        if load_geom:
            s = Scene()
            s.read(fgeom, 'BGEOM')
            geom = {sh.id: sh.geometry for sh in s}
            g.add_property('geometry')
            g.property('geometry').update(geom)

        return g

    def get_midribs(self, g, resample=False):
        visible_length = g.property('visible_length')
        blades = (vid for vid in g.vertices(scale=g.max_scale() - 1) if
                g.label(vid).startswith('blade'))
        vids = (vid for vid in blades if visible_length[vid] > 0)
        metamer = {vid: g.complex(vid) for vid in vids}
        axe = {vid: g.complex(metamer[vid]) for vid in metamer}
        plant = {vid: g.complex(axe[vid]) for vid in axe}
        p = g.property('species')
        species = {vid: p.get(plant[vid], 0) for vid in plant}

        midribs = {vid: self.leaves[species[vid]].midrib(g.node(vid), resample=resample) for
                   vid in species}
        #
        anchor = g.property('anchor_point')
        midribs_anchor = {
            vid: [anchor[cid] for cid in g.components(vid) if cid in anchor] for
            vid
            in midribs}
        hins = {k: v[0][2] + midribs[k][2] for k, v in
                midribs_anchor.items() if len(v) > 0}

        ntop = g.property('ntop')
        res = [pandas.DataFrame({'vid': vid,
                                 'ntop': ntop[vid],
                                 'metamer': int(
                                     g.label(metamer[vid]).split('metamer')[1]),
                                 'axe': g.label(axe[vid]),
                                 'plant': int(
                                     g.label(plant[vid]).split('plant')[1]),
                                 'species': species[vid],
                                 'x': midribs[vid][0],
                                 'y': midribs[vid][1],
                                 'hins': hins[vid]}) for vid in
               hins]  # hins keys are for midribs keys wich also have a geometry (anchor point

        return pandas.concat(res)

    def midrib_statistics(self, g):
        data = self.get_midribs(g)
        return midrib_statistics(data)
