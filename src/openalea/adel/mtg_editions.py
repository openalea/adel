# -*- python -*-
#
#       Copyright 2015 INRIA - CIRAD - INRA
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       WebSite : https://github.com/openalea-incubator/adel
#
# ==============================================================================
"""new mtg edition function (should be integrated in new mtg"""

from openalea.mtg.traversal import iter_mtg
from openalea.adel.newmtg import (
    internode_elements,
    sheath_elements,
    blade_elements,
    convert,
    properties_from_dict,
    adel_metamer,
)
from openalea.mtg import MTG, fat_mtg
from openalea.adel.exception import AdelDeprecationError


def find_label(label, g, constrained_in=0):
    """
    Find the vertex ids identified by label in g if it exists.
    Args:
        label: (str) the label to be matched
        g: a MTG
        constrained_in: a vid defining the complex the label is search


    Returns: (list of int) the vertex id of matching labels

    """
    ci = constrained_in
    labels = g.property("label")
    if ci == 0:
        return [k for k in g.vertices_iter() if labels.get(k, "") == label]
    else:
        # Hack car iter_mtg peut recycler sur le reste du mtg
        cscale = g.scale(ci)
        return [
            k
            for k in iter_mtg(g, ci)
            if (labels.get(k, "") == label) and (g.complex_at_scale(k, cscale) == ci)
        ]


def find_plants(g):
    """the vid of all plants in g

    Args:
        g: an mtg

    Returns: (list of int) the vid of plants in g

    """
    vids = g.component_roots_at_scale(g.root, 1)
    if len(vids) > 0:
        labels = g.property("label")
        return [vid for vid in vids if labels[vid].startswith("plant")]
    else:
        return vids


def find_metamers(g, plant_number=1, axe_label="MS"):
    """the vid of all metamers beared by axe indentified by its label

    Args:
        g: an mtg
        plant_number: (int) the number of the plant
        axe_label: (str) the label of the axe

    Returns:
        vid_plant (int) the vid of the plant bearing the metamers(int)
        vid_axe (int) the vid of the axe bearing the metamers
        vid_metamers (list of int) the vid of metamers

    """

    plant_label = "plant" + str(plant_number)
    vid_plant = find_label(plant_label, g)
    if len(vid_plant) > 0:
        vid_plant = vid_plant[0]
    else:
        raise ValueError("plant " + plant_label + " not found in g")

    plant_axes = g.components(vid_plant)
    vid_axe = [vid for vid in find_label(axe_label, g) if vid in plant_axes]
    if len(vid_axe) > 0:
        vid_axe = vid_axe[0]
    else:
        raise ValueError(
            "axis " + axe_label + " not found in g on plant " + plant_label
        )

    vids = g.components(vid_axe)
    if len(vids) > 0:
        labels = g.property("label")
        metamers = [vid for vid in vids if labels[vid].startswith("metamer")]
    else:
        metamers = []
    return vid_plant, vid_axe, metamers


def add_plant(
    g,
    plant_number=None,
    plant_properties=None,
    axis_properties=None,
    metamer_properties=None,
    collar_properties=None,
):
    """Add a plant identified by its number to an mtg.
    The plant is created with 5 sub-scales: plant / mainstem /
    metamer0 / collar / baseElement
    If the plant already exists, nothing is added

    Args:
        g: an mtg representing the canopy
        plant_number: (int) the number of the plant to add.
        If None a plant_number will be max(existing_plant_numbers) + 1
        plant_properties: a dict of properties associated to the plant vertex
        axis_properties: a dict of properties associated to the mainstem
        metamer_properties: a dict of properties associated to metamer0
        collar_properties:  a dict of properties associated to collar

    Returns: the vid of the plant created or found

    """
    plant_properties = plant_properties or {}
    axis_properties = axis_properties or {}
    metamer_properties = metamer_properties or {}
    collar_properties = collar_properties or {}

    if plant_number is None:
        labels = g.property("label")
        plant_labels = [
            labels[vid] for vid in g.component_roots_at_scale(g.root, scale=1)
        ]
        if len(plant_labels) > 0:
            plant_number = max([int(x.split("plant")[1]) for x in plant_labels]) + 1
        else:
            plant_number = 1

    label = "plant" + str(plant_number)
    found = find_label(label, g)
    if len(found) > 0:
        return found[0]
    else:
        vid_plant = g.add_component(
            g.root, label=label, edge_type="/", **plant_properties
        )
        vid_axe = g.add_component(
            vid_plant, edge_type="/", label="MS", **axis_properties
        )
        vid_metamer = g.add_component(
            vid_axe, edge_type="/", label="metamer0", **metamer_properties
        )
        vid_organ = g.add_component(
            vid_metamer, edge_type="/", label="collar", **collar_properties
        )
        vid_base = g.add_component(vid_organ, edge_type="/", label="baseElement")
        g.add_child(vid_base, edge_type="<", label="topElement")
        return vid_plant


def add_vegetative_metamer(
    g,
    plant_number=1,
    axe_label="MS",
    metamer_properties=None,
    internode_properties=None,
    sheath_properties=None,
    blade_properties=None,
):
    """Add a vegetative metatmer at the top of an axe

    Args:
        g: the MTG
        plant_number: (int) the plant number
        axe_label: the label of the axe bearing the new metamer
        metamer_properties: a dict of properties associated to the metamer
        internode_properties: a dict of properties associated to the internode
        sheath_properties: a dict of properties associated to the sheath
        blade_properties: a dict of properties associated to the blade

    Returns:
        the vid of the new metamer
    """

    metamer_properties = metamer_properties or {}
    internode_properties = internode_properties or {}
    sheath_properties = sheath_properties or {}
    blade_properties = blade_properties or {}

    vid_plants, vid_axe, metamers = find_metamers(g, plant_number, axe_label)
    num_metamer = len(metamers)  # start at zero
    label = "metamer" + str(num_metamer)
    # to check with christophe : I want the 'tip'
    vid_parent_metamer = max(g.components(vid_axe))
    vid_metamer = g.add_component(
        vid_axe, edge_type="/", label=label, **metamer_properties
    )
    g.add_child(vid_parent_metamer, child=vid_metamer, edge_type="<")
    # add organs
    # parent = internode or collar = first organ of the preceding metamer
    vid_parent_organ = min(g.components(vid_parent_metamer))
    # Top of the parent organ
    vid_parent_elt = find_label("topElement", g, vid_parent_organ)[0]
    # add internode
    new_organ = g.add_component(
        vid_metamer, label="internode", edge_type="/", **internode_properties
    )
    base_elt = g.add_component(new_organ, label="baseElement", edge_type="/")
    vid_parent_organ = g.add_child(vid_parent_organ, child=new_organ, edge_type="<")
    g.add_child(vid_parent_elt, child=base_elt, edge_type="<")
    vid_parent_elt = g.add_child(base_elt, label="topElement", edge_type="<")
    # sheath
    new_organ = g.add_component(
        vid_metamer, label="sheath", edge_type="/", **sheath_properties
    )
    base_elt = g.add_component(new_organ, label="baseElement", edge_type="/")
    vid_parent_organ = g.add_child(vid_parent_organ, child=new_organ, edge_type="+")
    g.add_child(vid_parent_elt, child=base_elt, edge_type="+")
    vid_parent_elt = g.add_child(base_elt, label="topElement", edge_type="<")
    # blade
    new_organ = g.add_component(
        vid_metamer, label="blade", edge_type="/", **blade_properties
    )
    base_elt = g.add_component(new_organ, label="baseElement", edge_type="/")
    g.add_child(vid_parent_organ, child=new_organ, edge_type="<")
    g.add_child(vid_parent_elt, child=base_elt, edge_type="<")
    g.add_child(base_elt, label="topElement", edge_type="<")

    return vid_metamer


def insert_elements(g, vid_organ, elements, before=None):
    """Insert elements between base and top element of an organ

    Args:
        g: the mtg
        vid_organ: the vertex id of the targetd organ
        elements: a list of dict with element properties
        before : the vertex id of the element before which elements should be
        inserted. If None (default) elements are inserted before top element

    Returns:
        a list of vid of the inserted elements
    """

    inserted = []
    if before is None:
        before = find_label("topElement", g, vid_organ)[0]
    for element in reversed(elements):
        before = g.insert_parent(before, edge_type="<", **element)
        inserted.append(before)
    return inserted


def add_axe(
    g,
    label,
    plant_number=1,
    axis_properties=None,
    metamer_properties=None,
    collar_properties=None,
):
    """Add an axe identified by its label to a plant identifed by its number.
    The axe is created with 4 sub-scales axe /metamer0 / collar / baseElement
    If the axe already exists, nothing is added
    If bearing metamers are missing they are added

    Args:
        g: an mtg representing the canopy
        label: (str) a string identifying the axe (eg 'T1' for an axe beared by
        metamer 1 of main stem, 'T0.1' for the axe beared by metamer1 of T0)
        plant_number: (int) the number of the plant bering the axe.
        axis_properties: a dict of properties associated to the mainstem
        metamer_properties: a dict of properties associated to metamer0
        collar_properties:  a dict of properties associated to collar

    Returns: the vid of the axe created or found

    """
    plant_label = "plant" + str(plant_number)
    vid_plant = find_label(plant_label, g)
    if len(vid_plant) > 0:
        vid_plant = vid_plant[0]
    else:
        raise ValueError("plant " + plant_label + " not found in g")

    axes = g.components(vid_plant)
    found = [vid for vid in find_label(label, g) if vid in axes]
    if len(found) > 0:
        return found[0]

    axis_properties = axis_properties or {}
    metamer_properties = metamer_properties or {}
    collar_properties = collar_properties or {}

    labels = g.property("label")
    axe_code = label.split("T")[1].split(".")

    if len(axe_code) > 1:
        parent_axe_label = "T" + ".".join(axe_code[:-1])
    else:
        parent_axe_label = "MS"

    # recursively add parent axes if needed
    found = [vid for vid in find_label(parent_axe_label, g) if vid in axes]
    if len(found) == 0:
        add_axe(g, parent_axe_label, plant_number)

    vid_axe = g.add_component(vid_plant, edge_type="/", label=label, **axis_properties)
    vid_metamer = g.add_component(
        vid_axe, edge_type="/", label="metamer0", **metamer_properties
    )
    vid_organ = g.add_component(
        vid_metamer, edge_type="/", label="collar", **collar_properties
    )
    vid_base_elt = g.add_component(vid_organ, edge_type="/", label="baseElement")
    vid_top_elt = g.add_component(vid_organ, edge_type="/", label="topElement")

    vid_plant, vid_parent_axe, parent_metamers = find_metamers(
        g, plant_number, parent_axe_label
    )
    bearing_metamer = int(axe_code[-1])
    last_metamer = len(parent_metamers) - 1
    if bearing_metamer > last_metamer:
        for i in range(last_metamer + 1, bearing_metamer + 1):
            new = add_vegetative_metamer(g, plant_number, labels[vid_parent_axe])
            parent_metamers.append(new)

    vid_parent_metamer = sorted(parent_metamers)[bearing_metamer]
    # parent organ = internode = first organ of the preceding metamer
    vid_parent_organ = min(g.components(vid_parent_metamer))
    # base of the parent organ = first element
    vid_parent_elt = find_label("baseElement", g, vid_parent_organ)[0]

    g.add_child(vid_parent_axe, child=vid_axe, edge_type="+")
    g.add_child(vid_parent_metamer, child=vid_metamer, edge_type="+")
    g.add_child(vid_parent_organ, child=vid_organ, edge_type="+")
    g.add_child(vid_parent_elt, child=vid_base_elt, edge_type="+")
    g.add_child(vid_base_elt, child=vid_top_elt, edge_type="<")

    return vid_axe


def update_organ_elements(g, leaves=None, split=False, phyllochron=None):
    """Set / update organ elements

    Args:
        g: an adel mtg
        leaves: a leaf shape database
        split: (bool) flag trigering the separation between senescent and green
        part of an organ

    Returns:

    """
    labels = g.property("label")
    length = g.property("length")
    visible_length = g.property("visible_length")
    rolled_length = g.property("rolled_length")
    senesced_length = g.property("senesced_length")
    azimuth = g.property("azimuth")
    inclination = g.property("inclination")
    diameter = g.property("diameter")
    sectors = g.property("n_sect")
    shape_mature_length = g.property("shape_mature_length")
    shape_max_width = g.property("shape_max_width")
    shape_key = g.property("shape_key")
    area = g.property("area")
    species = g.property("species")
    age = g.property("age")

    for organ in g.vertices(scale=4):
        if labels[organ].startswith("internode"):
            elts = internode_elements(
                length[organ],
                visible_length[organ],
                senesced_length[organ],
                azimuth[organ],
                inclination[organ],
                diameter[organ],
                split=split,
            )
        elif labels[organ].startswith("sheath"):
            elts = sheath_elements(
                length[organ],
                visible_length[organ],
                senesced_length[organ],
                azimuth[organ],
                inclination[organ],
                diameter[organ],
                split=split,
            )
        elif labels[organ].startswith("blade"):
            l = leaves[species[organ]]
            if l is not None and l.dynamic:
                lctype, lcindex, _ = shape_key[organ]
                axe = labels[g.complex(g.complex(organ))]
                age_index = l.get_age_index(
                    float(age[organ]) / phyllochron.get(axe, phyllochron["T1"]) - 0.3
                )
                shape_key[organ] = (lctype, lcindex, age_index)
            elts = blade_elements(
                sectors[organ],
                length[organ],
                visible_length[organ],
                rolled_length[organ],
                senesced_length[organ],
                shape_mature_length[organ],
                shape_max_width[organ],
                shape_key[organ],
                leaves=leaves[species[organ]],
                split=split,
            )
        else:
            elts = []

        if len(elts) > 0:
            # update area at organ scale
            if labels[organ].startswith("blade"):
                area[organ] = sum(
                    [
                        elt["area"]
                        for elt in elts
                        if elt["label"].startswith("LeafElement")
                    ]
                )
            else:
                area[organ] = sum(
                    [
                        elt["area"]
                        for elt in elts
                        if elt["label"].startswith("StemElement")
                    ]
                )

            # insert elts and updates at element scale
            if len(g.components(organ)) == 2:  # only top and base element
                insert_elements(g, organ, elts)
            else:
                # remove unmatched elt
                def remove_element(g, vid):
                    g.remove_vertex(vid, reparent_child=True)
                    # TODO hack waiting for bug correction in mtg
                    g._remove_vertex_properties(vid)

                newlabels = [e["label"] for e in elts] + ["baseElement", "topElement"]
                labels = g.property("label")
                for elt in g.components(organ):
                    if labels[elt] not in newlabels:
                        remove_element(g, elt)
                # insert new
                insertion_stack = []
                for elt in elts:
                    label = elt.pop("label")
                    found = find_label(label, g, organ)
                    if found:
                        vid_elt = found[0]
                        if len(insertion_stack) > 0:  # flush insertion stack
                            insert_elements(g, organ, insertion_stack, before=vid_elt)
                            insertion_stack = []
                        for k in elt:
                            g.property(k)[vid_elt] = elt[k]
                    else:
                        elt["label"] = label
                        insertion_stack.append(elt)
                if len(insertion_stack) > 0:  # no successor found
                    insert_elements(g, organ, insertion_stack)

    return g


def new_mtg_factory(
    parameters,
    metamer_factory=adel_metamer,
    leaf_sectors=1,
    leaves=None,
    stand=None,
    axis_dynamics=None,
    add_elongation=False,
    topology=("plant", "axe_id", "numphy"),
    split=False,
    aborting_tiller_reduction=1.0,
    leaf_db=None,
):
    """A 'clone' of mtg_factory that uses mtg_edition functions"""

    if leaf_db is not None:
        raise AdelDeprecationError(
            "leaf_db argument is deprecated, use leaves argument instead"
        )

    if leaves is None:
        dynamic_leaf_db = {0: False}
        leaves = {0: None}
    else:
        dynamic_leaf_db = {k: leaves[k].dynamic for k in leaves}

    g = MTG()

    # buffers
    # for detection of newplant/newaxe
    prev_plant = 0
    prev_axe = -1

    dp = parameters
    nrow = len(dp["plant"])
    species = 0
    for i in range(nrow):
        plant, num_metamer = [
            int(convert(dp.get(x)[i], undef=None))
            for x in [topology[e] for e in [0, 2]]
        ]
        axe = dp.get(topology[1])[i]
        args = properties_from_dict(dp, i, exclude=topology)

        # Add plant if new
        if plant != prev_plant:
            if axe != "MS":
                raise ValueError(
                    "Main stem is expected first when a new plant is declared"
                )

            position, azimuth = (0, 0, 0), 0
            if stand and len(stand) >= plant:
                position, azimuth = stand[plant - 1]
            plant_properties = {
                "position": position,
                "azimuth": azimuth,
                "refplant_id": args.get("refplant_id"),
                "species": species,
            }

            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            mainstem_properties = dict(
                timetable=timetable,
                HS_final=args.get("HS_final"),
                nff=args.get("nff"),
                hasEar=args.get("hasEar"),
                azimuth=args.get("az_insertion"),
            )

            add_plant(
                g,
                plant_number=plant,
                plant_properties=plant_properties,
                axis_properties=mainstem_properties,
            )
            prev_plant = plant
        # Add axis
        if axe != prev_axe and axe != "MS":
            timetable = None
            if axis_dynamics:
                timetable = axis_dynamics[str(plant)][str(axe)]
            axe_properties = dict(
                timetable=timetable,
                HS_final=args.get("HS_final"),
                nff=args.get("nff"),
                hasEar=args.get("hasEar"),
                azimuth=args.get("az_insertion"),
            )
            add_axe(g, axe, plant, axis_properties=axe_properties)
            prev_axe = axe
        # Add metamer
        assert num_metamer > 0
        # args are added to metamers only if metamer_factory is none,
        # otherwise compute metamer components
        components = []
        if metamer_factory:
            xysr_key = None
            if leaves[species] is not None and "LcType" in args and "LcIndex" in args:
                lctype = int(args["LcType"])
                lcindex = int(args["LcIndex"])
                if lctype != -999 and lcindex != -999:
                    age = None
                    if dynamic_leaf_db[species]:
                        age = (
                            float(args["rph"]) - 0.3
                        )  # age_db = HS - rank + 1 = ph - 1.3 - rank +1 = rph - .3
                        if age != "NA":
                            age = max(0, int(float(age)))
                    xysr_key = leaves[species].get_leaf_key(lctype, lcindex, age)

            elongation = None
            if add_elongation:
                startleaf = -0.4
                endleaf = 1.6
                stemleaf = 1.2
                startE = endleaf
                endE = startE + (endleaf - startleaf) / stemleaf
                endBlade = endleaf
                if args["Gl"] > 0:
                    endBlade = args["Ll"] / args["Gl"] * (endleaf - startleaf)
                elongation = {
                    "startleaf": startleaf,
                    "endBlade": endBlade,
                    "endleaf": endleaf,
                    "endE": endE,
                }
            if "ntop" not in args:
                args.update({"ntop": None})
            if "Gd" not in args:
                args.update({"Gd": 0.19})
            args.update({"split": split})
            if args.get("HS_final") < args.get("nff"):
                for what in (
                    "Ll",
                    "Lv",
                    "Lr",
                    "Lsen",
                    "L_shape",
                    "Lw_shape",
                    "Gl",
                    "Gv",
                    "Gsen",
                    "Gd",
                    "El",
                    "Ev",
                    "Esen",
                    "Ed",
                ):
                    args.update({what: args.get(what) * aborting_tiller_reduction})
            components = metamer_factory(
                Lsect=leaf_sectors,
                shape_key=xysr_key,
                elongation=elongation,
                leaves=leaves[species],
                **args,
            )
            args = {"L_shape": args.get("L_shape")}
        #
        metamer_properties = args
        internode, sheath, blade = None, None, None
        elts = {k: [] for k in ("internode", "sheath", "blade")}
        if len(components) > 0:
            internode, sheath, blade = components
            internode.pop("label")
            sheath.pop("label")
            blade.pop("label")
            elts["internode"] = internode.pop("elements")
            elts["sheath"] = sheath.pop("elements")
            elts["blade"] = blade.pop("elements")

        vid_metamer = add_vegetative_metamer(
            g,
            plant,
            axe,
            metamer_properties=metamer_properties,
            internode_properties=internode,
            sheath_properties=sheath,
            blade_properties=blade,
        )
        for organ in g.components(vid_metamer):
            label = g.property("label")[organ]
            if label in elts and len(elts[label]) > 0:
                insert_elements(g, organ, elts[label])
    return fat_mtg(g)
