def grow(n, time):
    """
    This function is called on each node to compute its evolution with time
    The following properties available as parameters:
        - on metamer: final_length of the leaf, start_tt of the leaf, end_tt of the leaf, start_insertion, end_insertion, insertion_angle, frate(falling rate in degrees.dd-1)
        - LeafElement: start_tt of the sector, end_tt of the sector

    The following properties are to be updated:
        - metamer: length of the leaf, start_tt, end_tt, insertion_angle
        - LeafElement: age, length( of the sector), color
        - StemElement : age,length,color

    """
    # Detect type of node (stem or leaf)
    if "Leaf" in n.label:
        # update metamer parameters (not to be changed)
        metamer = n.complex()
        if metamer.final_length is None:
            metamer.final_length = n.final_length
        metamer.start_insertion = 0
        metamer.end_insertion = n.Linc

        # compute age
        n.age = 0

        # colorise element
        n.color = 0, 180, 0

        # compute length of the entire leaf

        metamer.length = metamer.final_length

        # compute length of the sector
        n.length = (
            metamer.final_length
            * (n.end_tt - n.start_tt)
            / (metamer.end_tt - metamer.start_tt)
        )

        # Compute insertion angle of the leaf when passing on the first sector
        # dectect the first sector (not to be changed)
        if (n.start_tt <= time < n.end_tt) or (
            (time > metamer.end_tt) and n.edge_type() == "+"
        ):
            metamer.insertion_angle = metamer.end_insertion

    # fin  Leaf Element

    else:
        n.length = n.final_length
        n.age = 0
        n.color = 0, 90, 0
