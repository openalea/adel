from openalea.adel.mtg import *


def CanMTGInterpreter(symbols, lsystem_string):
    """Generates a CanMTG from geometric functions and an Lsystem string"""
    g = CanMTG(symbols, lsystem_string)
    # write the node code here.

    # return outputs
    return (g,)
