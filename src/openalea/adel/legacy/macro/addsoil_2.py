import string
import copy


def can_line(pts, label):
    return "p 1 %s 3 %s" % (
        str(label),
        " ".join("%.6f" % x for i in range(0, 3) for x in pts[i]),
    )


def addsoil_2(caribuscene, x0):
    """Add a soil along with the pattern of a CaribuScene (Supposing the pattern is aligned with coordinate axis)"""
    pat = caribuscene.pattern
    xy = list(map(string.split, pat.splitlines()))
    A = list(map(float, xy[0]))
    C = list(map(float, xy[1]))
    if A[0] > C[0]:
        A = list(map(float, xy[1]))
        C = list(map(float, xy[0]))
    if C[1] < A[1]:
        D = [A[0], C[1]]
        B = [C[0], A[1]]
    else:
        B = [A[0], C[1]]
        D = [C[0], A[1]]
    A.append(x0)
    B.append(x0)
    C.append(x0)
    D.append(x0)

    label = "000000000000"

    newscene = copy.copy(caribuscene)
    newscene.scene = "\n".join(
        [can_line((A, B, C), label), can_line((C, D, A), label), caribuscene.scene]
    )

    # return outputs
    return newscene
