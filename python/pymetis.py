import os
import copy
import datetime
import numpy as np
import six
from pydons import MatStruct
import importlib

_engine = None
_mlabbridge = None
_mpaths = []

try:
    from matlab.engine import start_matlab
    MATLAB_ENGINE = True
except ImportError:
    MATLAB_ENGINE = False
try:
    import pymatbridge
    PYMATBRIDGE = True
except ImportError:
    PYMATBRIDGE = False

if not (MATLAB_ENGINE or PYMATBRIDGE):
    raise ImportError('neither matlab nor pymatbridge available')


__all__ = ['metis', 'MetisX']


def engine():
    global _engine
    if _engine is None:
        _engine = start_matlab()
        _engine.addpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
        _engine.addpath(os.path.abspath(os.path.dirname(__file__)))
        for path in _mpaths:
            _engine.addpath(os.path.abspath(path))
        try:
            _engine.metis_path()
        except Exception:
            raise Exception('metis_path not found on MATLAB path - '
                            'use addpath or '
                            'include the path to METIS in your MATLABPATH')
    return _engine


def mlabbridge():
    global _mlabbridge
    if _mlabbridge is None:
        _mlabbridge = pymatbridge.Matlab(maxtime=120)
        _mlabbridge.start()
        _mlabbridge.run_func('addpath',
                     os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
        _mlabbridge.run_func('addpath', os.path.abspath(os.path.dirname(__file__)))
        for path in _mpaths:
            _mlabbridge.run_func('addpath', os.path.abspath(path))
        try:
            _mlabbridge.run_func('metis_path')
        except Exception:
            raise Exception('metis_path not found on MATLAB path - '
                            'use addpath or '
                            'include the path to METIS in your MATLABPATH')
    return _mlabbridge


def addpath(path):
    """Add path to MATLAB's path
    """
    global _mpaths
    _mpaths.append(os.path.abspath(path))
    # is the engine is already running we have to execute addpath
    if MATLAB_ENGINE and _engine is not None:
        _engine.addpath(os.path.abspath(path))
    elif PYMATBRIDGE and _mlabbridge is not None:
        _mlabbridge.run_func('addpath', os.path.abspath(path))


# TODO remote engine (via Pyro4?)


def metis(z0dinput, fast=True, convert=True, bridge=None, save_file=''):
    """Call METIS using Matlab engine

    Args:
        z0dinput (dict): z0dinput structure
        fast (bool): use fast calculation (zerodfast)
        convert (bool): convert result to numpy
        bridge (str): selects the matlab bridge: matlab or pymatbridge,
            None for automatic selection

    Returns:
        dict: post structure
    """

    if bridge is None:
        if MATLAB_ENGINE:
            bridge = 'matlab'
        elif PYMATBRIDGE:
            bridge = 'pymatbridge'

    if bridge.lower() == 'matlab':
        option = as_matlab(z0dinput.option)
        cons = as_matlab(z0dinput.cons)
        geo = as_matlab(z0dinput.geo)
        exp0d = as_matlab(z0dinput.exp0d)
        zerod, void_, profil0d = engine().zerod_wrap(option,
                                                     cons,
                                                     geo,
                                                     exp0d,
                                                     fast,
                                                     save_file,
                                                     nargout=3)

        post = MatStruct()
        post.zerod = MatStruct(zerod)
        post.profil0d = MatStruct(profil0d)
        post.z0dinput = MatStruct(z0dinput)

        if convert:
            post = mat2np(post)

        return post

    elif bridge.lower() == 'pymatbridge':
        zerod_wrap = os.path.join(os.path.dirname(__file__), 'zerod_wrap.m')
        # TODO (1, 1) ndarrays of strings are turned into cell arrays
        for k, v in z0dinput.option.items():
            if (isinstance(v, np.ndarray) and np.issubsctype(v, np.str_) and
                v.size == 1):
                z0dinput.option[k] = v.squeeze()[()]

        res = mlabbridge().run_func(zerod_wrap,
                                    z0dinput.option,
                                    z0dinput.cons,
                                    z0dinput.geo,
                                    z0dinput.exp0d,
                                    fast,
                                    save_file,
                                    nargout=3)

        post = MatStruct()
        post.zerod = MatStruct(res['result'][0])
        post.profil0d = MatStruct(res['result'][2])
        post.z0dinput = MatStruct(z0dinput)

        return post


class MetisX(object):
    """METIS as a function of a single input vector X.
    Particularly useful for optimizers.

    Example 1: fied time points, cons.ip input
        metisx = MetisX('file.mat', (0, (10, 20), 50, 100),
                        OrderedDict(('ip', ),
                                    ('pecrh', ))
        post = metisx([3e5, 1e6, 5e6])

    Args:
        input_data (str or dict): metis simulation file or z0dinput data structure
        x2metis_input_function (callable): function that converts an x-vector to the z0dinput Metis structure
        obj_func (callable): objective function getting the post Metis output as its argument
        obj_func_inp (callable): objective function calculated from inputs only (z0dinput)
            the value is added to cost_func so be careful not to include a single cost twice
        obj_func_limits (tuple): call metis only if obj_func_limits[0] < cost_func_inp(z0dinput) < cost_func_lim[1] (default no limit)
        memcache (bool): use memory cache (simple dict)
        cachedir (str): use disk cache (using joblib), None to disable
        matlab_bridge: user-supplied matlab bridge, None for automatic
    """

    def __init__(self,
                 input_data,
                 x2metis_input_function,
                 obj_func,
                 obj_func_inp=None,
                 obj_func_limits=(-np.inf, np.inf),
                 fast=True,
                 memcache=True,
                 cachedir=None,
                 matlab_bridge=None):
        super(MetisX, self).__init__()

        if isinstance(input_data, six.string_types):
            input_data = MatStruct.loadmat(input_data, 'post')
            self.z0dinput = mat2np(input_data.post.z0dinput)
        else:
            self.z0dinput = input_data

        # if cost functions are strings
        self.x2metis_input_function = import_function(x2metis_input_function)
        self.obj_func = import_function(obj_func)
        self.obj_func_inp = import_function(obj_func_inp)
        self.obj_func_limits = obj_func_limits
        self.fast = fast
        self.matlab_bridge = matlab_bridge
        self.memcache = memcache
        if memcache:
            self._cached_results = {}

        if cachedir is None:
            self.metis_call = self._x_metis_call
            self._joblibmem = None
        else:
            import joblib
            cachedir = os.path.join(cachedir,
                                    datetime.datetime.now().isoformat().replace(':', '-'))
            self._joblibmem = joblib.Memory(cachedir=cachedir, verbose=1)
            self.metis_call = self._joblibmem.cache(self._x_metis_call)

    # @staticmethod
    def _x_metis_call(self, x, kwargs, convert):
        # deepcopy prevents modifying the original input
        z0dinput = self.x2metis_input_function(copy.deepcopy(self.z0dinput), *x, **kwargs)
        post = metis(z0dinput, fast=self.fast, convert=convert, bridge=self.matlab_bridge)
        return post

    def __call__(self,
                 *x,
                 convert=True,
                 return_post=False,
                 **kwargs):

        post = None

        if self.memcache:
#             kwargs_hash = hash(frozenset(my_dict.items()))
#             x_hash = hash(x)
            inp_hash = (x, frozenset(kwargs.items()), self.fast, convert)
            post = self._cached_results.get(inp_hash)
            if post is not None:
                print('Result loaded from mem cache :-]')

        # if post is None:
        #     z0dinput = self.x2metis_input_function(self.z0dinput, *x, **kwargs)
        # else:
        #     z0dinput = deepcopy(post.z0dinput)

        # start objective function calculation
        obj = 0

#         # cost function of inputs
#         if self.obj_func_inp is not None:
#             obj = self.obj_func_inp(z0dinput)
#             # we have to keep the cost sign in the inequality
#             if obj > self.obj_func_limits[1]:
#                 if return_post:
#                     post = MatStruct()
#                     # TODO this must be a deepcopy
#                     post['z0dinput'] = deepcopy(z0dinput)
#                     post['error'] = (
#                         'obj_func_inp(self.z0dinput) > obj_func_inp, '
#                         'post cannot be calculated')
#                     return post
#                 else:
#                     return norm_cost(cost, self.ymin_norm, self.ymax_norm)

        # launch METIS
        if post is None:
            # post = self.metis_call(z0dinput, fast=fast, convert=convert, bridge=self.matlab_bridge)
            post = self.metis_call(x, kwargs, convert)

        if self.memcache and inp_hash not in self._cached_results:
            # TODO cache obj as well
            self._cached_results[inp_hash] = post
            print('Mem-cached the result, cache size: {}'.format(len(self._cached_results)))

        if return_post:
            res = post
        else:
            obj += self.obj_func(post)
            res = obj

        return res


def np2mat(x):
    """Convert NumPy arrays to Matlab variables.

    Args:
        x: data to be converted
    """
    import matlab

    if x.size == 0:
        return matlab.double([])
    x = np.asfortranarray(np.atleast_2d(x))
    # x = np.atleast_2d(x)
    xl = [g[()] for g in np.nditer(x, order='F')]
    is_complex = False
    if np.issubdtype(x.dtype, np.str_):
        # string arrays to plain str objects
        # TODO does not work for arrays of strings
        return str(x.squeeze())
    if np.can_cast(x.dtype, np.bool_):
        mtype = matlab.logical
    elif np.can_cast(x.dtype, np.int32):
        mtype = matlab.int32
    elif np.can_cast(x.dtype, np.int64):
        mtype = matlab.int64
    elif np.can_cast(x.dtype, np.float32):
        mtype = matlab.single
    elif np.can_cast(x.dtype, np.float64):
        mtype = matlab.double
    elif np.can_cast(x.dtype, np.complex128):
        mtype = matlab.double
        is_complex = True
    xm = mtype(xl, size=x.shape, is_complex=is_complex)
    return xm


def as_matlab(x):
    """Convert Pythod data to Matlab data.
    """
    import matlab

    if type(x) in (matlab.double, matlab.single, matlab.int8, matlab.int16,
                   matlab.int32, matlab.int64, matlab.uint8, matlab.uint16,
                   matlab.uint32, matlab.uint64, matlab.logical, matlab.object
                   ):
        return x
    elif isinstance(x, (np.ndarray, np.generic)):
        return np2mat(x)
    elif isinstance(x, int):
        # return matlab.int32(x)
        return x
    elif isinstance(x, float):
        # return matlab.double(x)
        return x
    elif isinstance(x, dict):
        return {k: as_matlab(v) for k, v in x.items()}
    elif isinstance(x, six.string_types):
        return x
    else:
        raise NotImplementedError('type {} not implemented'.format(type(x)))


def mat2np(x, squeeze=False):
    """Convert Matlab data to Python data.
    """
    if isinstance(x, dict):
        return MatStruct(((k, mat2np(v)) for k, v in x.items()))
    elif isinstance(x, (six.string_types, float, int, bool)):
        return x
    else:
        res = np.asfortranarray(x)
        if squeeze:
            res = res.squeeze()
        return res


def matfunc2py(matfunc, bridge=None, nargout=1):
    """Converts matlab function to Python function.
    """

    if bridge is None:
        if MATLAB_ENGINE:
            bridge = 'matlab'
        elif PYMATBRIDGE:
            bridge = 'pymatbridge'

    def pyfunc(*args, **kwargs):
        """Python wrapper for Matlab {} function
        """.format(matfunc)

        if bridge.lower() == 'matlab':
            args_mat = [as_matlab(arg) for arg in args]
            for key, value in kwargs.items():
                args_mat.append(str(key))
                args_mat.append(as_matlab(value))

            res_mat = getattr(engine(), matfunc)(*args_mat, nargout=nargout)
            if nargout == 1:
                res_py = mat2np(res_mat)
            else:
                res_py = (mat2np(r) for r in res_mat)

            return res_py

        elif bridge.lower() == 'pymatbridge':
            wh = mlabbridge().run_code("which('{}')".format(matfunc))
            if (not wh['success']) or ('not found' in wh['content']['stdout']):
                raise Exception('Matlab function {} not found'.format(matfunc))
            matfunc_path = wh['content']['stdout'].strip()

            res = mlabbridge().run_func(matfunc_path, *args, nargout=nargout)

            return res

    return pyfunc


def import_function(func):
    """Import a function possible specified by a string
    """

    if isinstance(func, six.string_types):
        module_name = '.'.join(func.split('.')[:-1])
        func_name = func.split('.')[-1]
        if module_name:
            module = importlib.import_module(module_name)
            func = getattr(module, func_name)
        else:
            func = eval(func_name)

    return func


def apply_waveform(time_base, input_data, times, values):
    """Applies a linerly interpolated waveform
    """
    if times[0] > time_base[0]:
        times = np.hstack((time_base[0], times))
        values = np.hstack((input_data[0], values))
    if times[-1] < time_base[-1]:
        times = np.hstack((times, time_base[-1]))
        values = np.hstack((values, input_data[-1]))
    if np.iscomplexobj(values):
        res = np.atleast_2d(np.interp(time_base, times, values.real) +
                            1j * np.interp(time_base, times, values.imag))
    else:
        res = np.atleast_2d(np.interp(time_base, times, values))
    return res
