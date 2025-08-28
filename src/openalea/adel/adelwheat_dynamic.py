"""Prototype adel model that uses mtg edition functions"""

from openalea.mtg import MTG
from openalea.adel.astk_interface import AdelWheat
from openalea.adel.mtg_editions import (
    find_metamers,
    add_plant,
    add_vegetative_metamer,
    new_mtg_factory,
    update_organ_elements,
)
from openalea.adel.AdelR import plantSample
from openalea.adel.mtg_interpreter import mtg_interpreter, transform_geom


class AdelWheatDyn(AdelWheat):
    def build_stand(self):
        g = MTG()
        plants = self.axeT()
        sample = [int(p) for p in plantSample(self.pars)]
        for i, plant in plants.groupby("plant"):
            plant_properties = {
                "position": self.positions[i - 1],
                "azimuth": self.plant_azimuths[i - 1],
                "refplant_id": sample[i - 1],
            }
            ms = plant.loc[plant["axe"] == "MS", :].to_dict("list")
            ms_properties = {
                "HS_final": float(ms["HS_final"][0]),
                "nff": int(ms["nf"][0]),
                "hasEar": bool(int(ms["hasEar"][0])),
                "azimuth": float(ms["azTb"][0]),
            }
            add_plant(
                g, i, plant_properties=plant_properties, axis_properties=ms_properties
            )
        return g

    def add_metamer(self, g, plant=1, axe="MS"):
        vid_plant, vid_axe, metamers = find_metamers(g, plant, axe)
        nff = g.property("nff")[vid_axe]
        num_metamer = len(metamers)
        df = self.phytoT(axe)
        m = df.loc[df["n"] == num_metamer, :].to_dict("list")
        metamer_properties = {"L_shape": m["Ll"][0]}
        ntop = nff - num_metamer + 1
        shape_key = None
        lctype = int(m["Lindex"][0])
        lcindex = int(m["Lseed"][0])
        if lctype != -999 and lcindex != -999:
            shape_key = self.leaves[0].get_leaf_key(lctype, lcindex, age=None)

        internode_properties = {
            "ntop": ntop,
            "length": 0,
            "visible_length": 0,
            "senesced_length": 0,
            "diameter": m["Ed"][0],
            "azimuth": m["Azim"][0],
            "inclination": 0,
        }
        sheath_properties = {
            "ntop": ntop,
            "length": 0,
            "visible_length": 0,
            "senesced_length": 0,
            "diameter": m["Gd"][0],
            "azimuth": 0,
            "inclination": 0,
        }
        blade_properties = {
            "ntop": ntop,
            "length": 0,
            "visible_length": 0,
            "rolled_length": 0,
            "senesced_length": 0,
            "n_sect": self.nsect,
            "shape_mature_length": m["Ll"][0],
            "shape_max_width": m["Lw"][0],
            "shape_key": shape_key,
            "inclination": 0,
            "species": 0,
        }

        return add_vegetative_metamer(
            g,
            plant,
            axe,
            metamer_properties,
            internode_properties,
            sheath_properties,
            blade_properties,
        )

    def build_mtg(self, parameters, stand, **kwds):
        """temporary overwrite adel default"""
        g = new_mtg_factory(
            parameters,
            stand=stand,
            leaf_sectors=self.nsect,
            leaves=self.leaves,
            split=self.split,
            **kwds,
        )
        g = mtg_interpreter(g, self.leaves, min_length=self.min_length, classic=self.classic, face_up=self.face_up)
        return g

    def update_geometry(
        self, g, SI_units=False, properties_to_convert={"lengths": [], "areas": []}
    ):
        """Update MTG geometry.

        :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update the geometry.
            - `SI_units` (:class:`bool`) - A boolean indicating whether the MTG properties are expressed in SI units.
            - `properties_to_convert` (:class:`dict` of :class:`pandas.DataFrame`) - A dictionnary with the list of length properties area properties to be converted.
        :Returns:
            MTG with updated geometry
        :Returns Type:
            :class:`openalea.mtg.mtg.MTG`
        """

        if SI_units:
            self.convert_to_ADEL_units(g, properties_to_convert)

        # update elements
        g = update_organ_elements(g, self.leaves, self.split, self.phyllochron())
        g = mtg_interpreter(g, self.leaves, min_length=self.min_length, face_up=self.face_up, classic=self.classic)
        pos = g.property("position")
        az = g.property("azimuth")
        geom = g.property("geometry")
        for i, vid in enumerate(g.vertices(1)):
            pos[vid] = self.positions[i]
            az[vid] = self.plant_azimuths[i]
            for gid in g.components_at_scale(vid, g.max_scale()):
                if gid in geom:
                    geom[gid] = transform_geom(
                        geom[gid], self.positions[i], self.plant_azimuths[i]
                    )
        return g

    def convert_to_ADEL_units(self, g, properties_to_convert):
        """Converts the MTG to ADEL units. From m to cm for length properties and from m2 to cm2 for area properties.

        :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update.
            - `properties_to_convert` (:class:`dict` of :class:`pandas.DataFrame`) - A dictionnary with the list of length properties area properties to be converted.
        """
        for length_property in properties_to_convert["lengths"]:
            for vid, length in g.properties()[length_property].iteritems():
                g.properties()[length_property][vid] = length * 100  # Converts m to cm

        for area_property in properties_to_convert["areas"]:
            for vid, area in g.properties()[area_property].iteritems():
                g.properties()[area_property][vid] = area * 10000  # Converts m2 to cm2

    def convert_to_SI_units(self, g, properties_to_convert):
        """Converts the MTG to SI units from ADEL. From cm to m for length properties and from cm2 to m2 for area properties.

        :Parameters:
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update.
            - `properties_to_convert` (:class:`dict` of :class:`pandas.DataFrame`) - A dictionnary with the list of length properties area properties to be converted.
        """
        for length_property in properties_to_convert["lengths"]:
            for vid, length in g.properties()[length_property].iteritems():
                g.properties()[length_property][vid] = length / 100  # Converts cm to m

        for area_property in properties_to_convert["areas"]:
            for vid, area in g.properties()[area_property].iteritems():
                g.properties()[area_property][vid] = area / 10000  # Converts cm2 to m2
