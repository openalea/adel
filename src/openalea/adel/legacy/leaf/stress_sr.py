from .curvature import curvilinear_abscisse
import numpy as np


class Interval:
    def __init__(self, min_, max_):
        self._min = min_
        self._max = max_

    @property
    def min(self):
        return self._min

    @property
    def max(self):
        return self._max

    def __contains__(self, v):
        return self.min <= v <= self.max

    def __hash__(self):
        return hash((self.min, self.max))

    def __eq__(self, obj):
        return hash(obj) == hash((self.min, self.max))

    def __cmp__(self, obj):
        # incomplete implementation
        if isinstance(obj, Interval):
            if self.max <= obj.min:
                return -1
            elif self.min >= obj.max:
                return 1
            elif self == obj:
                return 0

    def __repr__(self):
        return "Interval(" + str(self.min) + "," + str(self.max) + ")"


def sort_intervals(d):
    intervals = dict(list(zip((Interval(m, M) for m, M in d), iter(d.values()))))

    keys = list(sorted(intervals))
    # add missing intervals with factor 1
    i = keys[0]
    if i.min > 0:
        intervals[Interval(0, i.min)] = 1
    for index in range(len(keys) - 1):
        i, j = keys[index : index + 2]
        assert i < j
        if i.max < j.min:
            intervals[Interval(i.max, j.min)] = 1
    j = keys[-1]
    if j.max < 1.0:
        intervals[Interval(j.max, 1)] = 1

    keys = list(sorted(intervals))
    return keys, intervals


def stress_sr(s, r, d):
    """ """
    keys, intervals = sort_intervals(d)

    i = 0
    interval = keys[i]
    factor = intervals[interval]
    p_min = interval.min
    # p_min
    l = []
    for p in s:
        if p not in interval:
            # accumulate values by taking account of the size of the previous interval modified by the factor
            p_min += (interval.max - interval.min) * factor
            i += 1
            interval = keys[i]
            factor = intervals[interval]
        l.append(p_min + (p - interval.min) * factor)

    new_s = np.array(l)

    return (
        new_s,
        r,
    )


# a changer avec passage par curvature ?
def stress_xy(s, x, y, d):
    """ """
    # import numpy as np
    from scipy.stats._support import unique

    keys, intervals = sort_intervals(d)

    s_new = [keys[0].min] + [k.max for k in keys]
    s_new = list(sorted(set(s).union(s_new)))
    x_new = np.interp(s_new, s, x)
    y_new = np.interp(s_new, s, y)

    i = 0
    interval = keys[i]
    factor = intervals[interval]
    p_min = interval.min
    x_min = x_new[0]
    y_min = y_new[0]

    dx = np.diff(x_new)
    dy = np.diff(y_new)

    # p_min
    l = [p_min]
    xl = []
    yl = []
    for index, p in enumerate(s[1:]):
        if p not in interval:
            # accumulate values by taking account of the size of the previous interval modified by the factor
            p_min += (interval.max - interval.min) * factor
            i += 1
            interval = keys[i]
            factor = intervals[interval]

        l.append(p_min + (p - interval.min) * factor)

        xl.append(dx[index] * factor)
        yl.append(dy[index] * factor)

    new_x = np.cumsum([x_min] + xl)
    new_y = np.cumsum([y_min] + yl)
    new_s = np.array(l)

    return new_s, new_x, new_y


def stress_xysr(x, y, s, r, d):
    """ """
    # import numpy as np
    from scipy.stats._support import unique

    keys, intervals = sort_intervals(d)

    new_s, r = stress_sr(s, r, d)

    s_xy = curvilinear_abscisse(x, y)
    s_max = s_xy.max()
    s_xy /= s_max
    new_x = x / s_max
    new_y = y / s_max

    s_xy, new_x, new_y = stress_xy(s_xy, new_x, new_y, d)

    xy = np.array([x, y])
    xy = unique(xy.T).T
    new_x, new_y = xy[0], xy[1]
    return new_x, new_y, new_s, r
