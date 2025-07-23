def AdelRunOptions(
    leafDuration,
    fracLeaf,
    stemDuration,
    dHS_col,
    dHS_en,
    senescence_leaf_shrink,
    epsillon,
    HSstart_inclination_tiller,
    rate_inclination_tiller,
    drop_empty,
):
    """Node interface for Adel General Option dict
    -leafDuration : the duration of leaf extension (blade and sheath); default value =2 phyllochron
    fracLeaf : the ratio between the duration of leaf extension before emergence (hidden part) and the duration of total extension of the leaf  ; default value =(0.4 )/(0.4+1.6)=0.2


    -fracLeaf =  leafDuration(hidden part) /leafDuration (total)


    -stemDuration: the duration of internode extension ; default value =1.66 phyllochron
    stemDuration= leafDuration/1.2


    -dHS_col: the phyllochronic delay between the time when {HS=n} and {the time of collar appearance of leaf n}; dHS_col is used in Adel to compute the HS of tiller required to determine the tiller inclination kinetics.
    (e.g.when collar of leaf n appear, the HS= n - dHS_col)


    -dHS_en: the phyllochronic delay between the end of leaf extension (blade and sheath) and the start of internode elongation  ; default value =0 phyllochron


    -senescence_leaf_shrink : the reduction factor applied to blade width of senescent part


    -epsillon: the minimum length taking into account in adel ; default value =1e-6 (cm)
    (this parameter was introduced to debug Caribu that didn't support small object )


    -HSstart_inclination_tiller': the haun stage of tiller at witch the stem start their inclination; default value =1
    {before this stage , the inclination angle is taken as constant (3 degree)}


    -rate_inclination_tiller: the rate of tiller inclination per HS; default value =30 degree //if (HS_tiller  =< HSstart_inclination_tiller) ,  TillerInclinationAngle=3 degree; ifelse  TillerInclinationAngle=Min (incT,rate_inclination_tiller*(HS_tiller- HSstart_inclination_tiller)  )


    -drop_empty: When checked allow to get an mtg with only the existing object (not empty object) at a precise time step.
    """
    adel_general_option_dict = {
        "leafDuration": leafDuration,
        "fracLeaf": fracLeaf,
        "stemDuration": stemDuration,
        "dHS_col": dHS_col,
        "dHS_en": dHS_en,
        "senescence_leaf_shrink": senescence_leaf_shrink,
        "epsillon": 1e-6,
        "HSstart_inclination_tiller": HSstart_inclination_tiller,
        "rate_inclination_tiller": rate_inclination_tiller,
        "drop_empty": drop_empty,
    }
    return (adel_general_option_dict,)
