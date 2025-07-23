try:
    from scipy import stats
except:
    from scipy import stats
import numpy as np
import pandas


def pandadf2adeldict(df):
    """convertit un dataframe panda en dictionaire de vecteur numpy"""
    d = df.to_dict()
    return dict((k, np.array([v for v in dv.values()])) for k, dv in d.items())


def dimension_fitting(id_axe_ref, csvdata_to_fit):
    def fit_length(tt, L_organ, prim_group):
        """
        tt: tip value (float)
        L_organ: the type of organ (str)
        prim_group: the primary group to use for fitting (DataFrame)
        return: the fitted length for tt (float)
        """

        for i in range(prim_group.index.size):
            max_dtt = prim_group["dtt_em_phytomer"][i]
            min_dtt = prim_group["dtt_em_phytomer"][i - 1]
            max_L = prim_group[L_organ][i]
            min_L = prim_group[L_organ][i - 1]

            if max_dtt < tt:
                pass
            elif max_dtt == tt:
                return prim_group[L_organ][i]
            elif max_dtt > tt:
                return ((max_L - min_L) * (tt - min_dtt) / (max_dtt - min_dtt)) + min_L

        return ((max_L - min_L) / (max_dtt - min_dtt)) * tt + (
            (max_L * min_dtt) - (min_L * max_dtt)
        ) / (min_dtt - max_dtt)

    def fit_width(tt, slope, intercept):
        """
        tt :tip value (float)
        slope: rate (float)
        return: fitted width (float)
        """
        return slope * tt + intercept

    def complete_df(df, id_dim_p):
        """
        df: table to fit (dataframe)
        id_dim_p: index of the primary axis used for the fit (int)
        return: table fitted (dataframe)
        """

        def select_from_id_dim(i):
            if df["id_dim"][i] == id_dim:
                return True
            return False

        id_dim = id_dim_p
        prim_group = df.select(select_from_id_dim)

        all_groups = df.groupby("id_dim").groups

        fitted_df = pandas.DataFrame()

        for id_dim_curr, line_numbers in all_groups.items():
            id_dim = id_dim_curr
            df_to_fit = df.select(select_from_id_dim)
            # select secondary axis: we supose that primary axis have an id_dim =(100+FinalLeafNumber)
            if id_dim_curr >= 200:
                for column_name in df.columns:
                    # fit the lengths
                    if column_name.startswith("L_"):
                        for i in range(df_to_fit.index.size):
                            tt = df_to_fit["dtt_em_phytomer"][i]
                            df_to_fit[column_name][i] = fit_length(
                                tt, column_name, prim_group
                            )
                    # fit the width
                    elif column_name.startswith("W_"):
                        first_W = prim_group[column_name][0]
                        last_W = prim_group[column_name][-1]
                        first_dtt = df_to_fit["dtt_em_phytomer"][0]
                        last_dtt = df_to_fit["dtt_em_phytomer"][-1]
                        # utiliser la fonction de regression lineaire pour determiner la pente et l'ordonne a l.origine
                        slope, intercept, r_value, p_value, std_err = stats.linregress(
                            [first_dtt, last_dtt], [first_W, last_W]
                        )

                        for i in range(df_to_fit.index.size):
                            tt = df_to_fit["dtt_em_phytomer"][i]
                            df_to_fit[column_name][i] = fit_width(tt, slope, intercept)

            fitted_df = fitted_df.append(df_to_fit)

        return fitted_df

    dimT_df = pandas.read_csv(csvdata_to_fit)
    dataframe_fitted = complete_df(dimT_df, id_axe_ref)
    del dataframe_fitted["dtt_em_phytomer"]
    #####################################Choisir le format de sortie#############################################

    #### output est un csv avec csvdata_fitted est le chemain d.accee plus le nom de fichier csv de sortie. exp(csvdata_fitted:D/dataframe_fitted.csv)
    # fitted_df.to_csv(csvdata_fitted, na_rep='NA', index=False)

    #### output est un pandasDataFrame
    # dataframe_fitted

    ####output est un np.ndarray:transforme pandas.dataframe en numpy array
    # nparray_fitted =dataframe_fitted.as_matrix(columns=(['id_dim','index_phytomer','L_blade','W_blade','L_sheath','W_sheath','L_internode','W_internode']))

    ####output est un dictionnaire:transforme pandas.dataframe en dictionnaire
    # dict_fitted=dataframe_fitted.to_dict()

    return ({"dimT": pandadf2adeldict(dataframe_fitted)},)
