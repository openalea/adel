import base64
import json
import numpy


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        if input object is a ndarray it will be converted into a dict holding dtype, shape and the data base64 encoded
        """
        if isinstance(obj, numpy.ndarray):
            data_b64 = base64.b64encode(obj.data)
            return dict(__ndarray__=data_b64, dtype=str(obj.dtype), shape=obj.shape)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder(self, obj)


def json_numpy_obj_hook(dct):
    """
    Decodes a previously encoded numpy ndarray
    with proper shape and dtype
    :param dct: (dict) json encoded ndarray
    :return: (ndarray) if input was an encoded ndarray
    """
    if isinstance(dct, dict) and "__ndarray__" in dct:
        data = base64.b64decode(dct["__ndarray__"])
        return numpy.frombuffer(data, dct["dtype"]).reshape(dct["shape"])
    return dct


# Overload dump/load to default use this behavior.
def dumps(*args, **kwargs):
    kwargs.setdefault("cls", NumpyEncoder)
    return json.dumps(*args, **kwargs)


def loads(*args, **kwargs):
    kwargs.setdefault("object_hook", json_numpy_obj_hook)
    return json.loads(*args, **kwargs)


def dump(*args, **kwargs):
    kwargs.setdefault("cls", NumpyEncoder)
    return json.dump(*args, **kwargs)


def load(*args, **kwargs):
    kwargs.setdefault("object_hook", json_numpy_obj_hook)
    return json.load(*args, **kwargs)
