from rpy import *
from numpy import transpose

data = r.load(".RData")

nerva = r.get("nerva")
nervj = r.get("nervj")
SRa = r.get("SRa")
SRj = r.get("SRj")

nervj["CPJ"] = nervj["CPJ010"]
del nerva["menuet"]

genotype = "F36"
gen = genotype

# Extract geometric information for each leaf
# database structure:
# genotype, plant number, leaf rank, data
# database key: rank, value : x,y,s,r


def build_db(gen, nervj, nerva, SRj, SRa):
    plant_id = set(nervj[gen].keys())
    plant_id.intersection_update(list(nerva[gen].keys()))
    plant_id.intersection_update(list(SRj[gen].keys()))
    plant_id.intersection_update(list(SRa[gen].keys()))
    leaves = {}
    for pid in plant_id:
        nerv = nervj[gen][pid]
        rad = SRj[gen][pid]
        rank_id = set(nerv.keys())
        rank_id.intersection_update(list(rad.keys()))
        for k in rank_id:
            x, y = transpose(nerv[k])
            s, r = transpose(rad[k])
            leaves.setdefault(k, []).append((x, y, s, r))

        nerv = nerva[gen][pid]
        rad = SRa[gen][pid]
        rank_id = set(nerv.keys())
        rank_id.intersection_update(list(rad.keys()))
        for k in rank_id:
            x, y = transpose(nerv[k])
            s, r = transpose(rad[k])
            leaves.setdefault(k, []).append((x, y, s, r))

    return leaves


if __name__ == "__main__":
    leaves = build_db(gen, nervj, nerva, SRj, SRa)
    import pickle as Pickle

    f = open("leaves.db", "w")
    Pickle.dump(leaves, f)
    f.close()
