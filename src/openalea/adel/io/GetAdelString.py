import os


def GetAdelString(lsysdir, iternumber):
    """Reurn the content of a numbered lsystem string file"""

    fn = "AleaStr%s.txt" % (iternumber)
    ff = open(lsysdir + "/" + fn)
    res = ff.read()
    ff.close()

    return res
