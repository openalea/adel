def grow(n, time):
    """
    This function is called on each node to compute its evolution with time
    The following properties available as parameters:
        - on metamer: final_length of the leaf, start_tt of the leaf, end_tt of the leaf, start_insertion, end_insertion, insertion_angle, frate(falling rate in degrees.dd-1)
        - LeafElement: start_tt of the sector, end_tt of the sector

    The following properties are to be updated:
        - metamer: length of the leaf, start_tt, end_tt, insertion_angle
        - LeafElement: age, length of the sector
        - StemElement : age,length

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
        n.age = time - n.start_tt

        # colorise as a function of age
        if n.age >= 350:
            if time - n.end_tt >= 350:
                n.color = 255, 0, 0
            else:
                n.color = 255, 170, 0
        else:
            n.color = 0, 180, 0

        # compute length of the entire leaf
        development_stage = (min(time, metamer.end_tt) - metamer.start_tt) / (
            metamer.end_tt - metamer.start_tt
        )
        metamer.length = metamer.final_length * development_stage

        # compute length of the sector
        n.length = (
            metamer.final_length
            * (min(time, n.end_tt) - n.start_tt)
            / (metamer.end_tt - metamer.start_tt)
        )

        # Compute insertion angle of the leaf when passing on the first sector
        # detect the first sector (not to be changed)
        if (n.start_tt <= time < n.end_tt) or (
            (time > metamer.end_tt) and n.edge_type() == "+"
        ):
            if time <= metamer.end_tt:
                metamer.insertion_angle = (
                    metamer.start_insertion
                    + development_stage
                    * (metamer.end_insertion - metamer.start_insertion)
                )
            else:
                metamer.insertion_angle = min(
                    90, metamer.end_insertion + metamer.frate * (time - metamer.end_tt)
                )

        # fin  Leaf Element

    else:
        n.length = (
            n.final_length
            * (min(time, n.end_tt) - n.start_tt)
            / (n.end_tt - n.start_tt)
        )
        n.age = time - n.start_tt
        n.color = 0, 90, 0
