from openalea.adel.AdelR import R_xydb, R_srdb


def extract_leaf_info(rdata_xy, rdata_sr):
    """
    Extract geometric information for each leaf database structure:
    genotype, plant number, leaf rank, data
    database key: rank, value : x,y,s,r
    """

    xy = R_xydb(rdata_xy)
    sr = R_srdb(rdata_sr)

    rank = list(xy.keys())
    rank.sort()
    # Assert that the two databases are compatible.
    rank_bis = list(sr.keys())
    rank_bis.sort()
    assert rank == rank_bis

    for k in rank:
        for i in range(len(xy[k])):
            xy[k][i] = (xy[k][i][0], xy[k][i][1], sr[k][0], sr[k][1])
    return xy


def test():
    rdata_xy = "data/So99.RData"
    rdata_sr = "data/SRSo.RData"
    leaves = extract_leaf_info(rdata_xy, rdata_sr)
    rank = list(leaves.keys())
    rank.sort()
    assert rank == ["1", "2", "3", "4", "5", "6"]
    return leaves


if __name__ == "__main__":
    leaves = test()
    import pickle as Pickle

    f = open("leaves_wheat.db", "w")
    Pickle.dump(leaves, f)
    f.close()
