def AND(xin, bool):
    bool = [x - 2 for x in bool]
    bool = [abs(x) for x in bool]
    xout = [x * y for x, y in zip(bool, xin)]
    return (xout,)
