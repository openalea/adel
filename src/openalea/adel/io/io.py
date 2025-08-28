# Standard import

import pickle as Pickle

from rpy2 import robjects

r = robjects.r

# needed for wralea


def load_leaf_data(fn):
    """Load leaf data obtained by measurement."""
    try:
        f = open(fn)
        leaves = Pickle.load(f)
    finally:
        f.close()

    return (leaves,)
