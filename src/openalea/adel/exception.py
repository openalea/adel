"""Defines Exception classes for adel"""


class AdelError(Exception):
    pass


class AdelParameterisationError(AdelError):
    pass


class AdelDeprecationError(AdelError):
    pass


class AdelImplementationError(AdelError):
    pass
