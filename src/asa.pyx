# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Adaptive Simulated Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@gmail.com>

"""Adaptive Simulated Annealing"""

import sys

import numpy as np
cimport numpy as np

from cpython cimport Py_INCREF, Py_DECREF

from asa_def cimport *

np.import_array()


asa_codes = {

    NORMAL_EXIT: """Given the criteria set, the search has run its
    normal course, e.g., when Limit_Acceptances or Limit_Generated is
    reached.""", 

    P_TEMP_TOO_SMALL: """A parameter temperature was too small using the
    set criteria. Often this is an acceptable status code.""",

    C_TEMP_TOO_SMALL: """The cost temperature was too small using the
    set criteria. Often this is an acceptable status code.""",

    COST_REPEATING: """The cost function value repeated a number of
    times using the set criteria. Often this is an acceptable status
    code.""",

    TOO_MANY_INVALID_STATES: """Too many repetitive generated states
    were invalid using the set criteria. A result of raising
    CostParameterError to include constraints.""",

    IMMEDIATE_EXIT: """The user has set OPTIONS−>Immediate_Exit to TRUE,
    or deleted file asa_exit_anytime when ASA_EXIT_ANYTIME is TRUE. This
    should not happen as it is converted to an exception.""",

    INVALID_USER_INPUT: """The user has introduced invalid input. When
    entering asa(), a function asa_test_asa_options() checks out many
    user−defined parameters and OPTIONS, and prints out invalid OPTIONS
    when ASA_PRINT is set to TRUE.""",

    INVALID_COST_FUNCTION: """The user has returned a value of the cost
    function to asa () which is not a valid number, e.g., not between
    -MAX_DOUBLE and MAX_DOUBLE.  Or, the user has returned a value of a
    parameter no longer within its proper range (excluding cases where
    the user has set the lower bound equal to the upper bound to remove
    a parameter from consideration).""",

    INVALID_COST_FUNCTION_DERIV: """While calculating numerical cost
    derivatives, a value of the cost function was returned which is not
    a valid number, e.g., not between -MAX_DOUBLE and MAX_DOUBLE. Or,
    while calculating numerical cost derivatives, a value of a parameter
    no longer within its proper range (excluding cases where the user
    has set the lower bound equal to the upper bound to remove a
    parameter from consideration) was set.""",

    CALLOC_FAILED: """Calloc memory allocation has failed in asa.c. This
    error will call Exit_ASA(), the location will be printed to stdout,
    and asa() will return the double -1 to the calling program. In
    asa_usr.c, if asa() returns this *exit_code a warning will be
    printed both to stdout and to USER_OUT.  Note that if a calloc
    memory allocation fails in asa_usr.c, this error will call
    Exit_USER() to print the location to stdout and then return -2.""",

}


class CostParameterError(Exception):
    """Exception to be raised in the cost function 
    if the parameters supplied are invalid."""
    pass


def asa(object func not None,
        np.ndarray[np.double_t, ndim=1] x0 not None,
        np.ndarray[np.double_t, ndim=1] xmin not None,
        np.ndarray[np.double_t, ndim=1] xmax not None,
        int full_output=False,
        tuple args=(),
        dict kwargs={},
        np.ndarray[np.int_t, ndim=1] parameter_type=None,
        long rand_seed=696969,
        int limit_acceptances=1000,
        int limit_generated=99999,
        int limit_invalid_generated_states=1000,
        double accepted_to_generated_ratio=1e-4,
        double cost_precision=1e-18,
        int maximum_cost_repeat=5,
        int number_cost_samples=5,
        double temperature_ratio_scale=1e-5,
        double cost_parameter_scale_ratio=1.,
        double temperature_anneal_scale=100., 
        int include_integer_parameters=False,
        int user_initial_parameters=False,
        int sequential_parameters=-1,
        double initial_parameter_temperature=1.,
        int acceptance_frequency_modulus=100,
        int generated_frequency_modulus=10000,
        int reanneal_cost=1,
        int reanneal_parameters=1,
        double delta_x=1e-3, 
        #asa_out_file="asa.log",
        ):
    """Adaptive Simulated Annealing

    Adaptive Simulated Annealing (ASA) is a C-language code developed to
    statistically find the best global fit of a nonlinear constrained
    non-convex cost-function over a D-dimensional space [1]_.
    
    Parameters
    -----------
    func : function
        Cost/merit function to minimize. Call signature:
            func(x, *args, **kwargs)
        Must return a float.
        Should raise CostParameterError if the parameters passed
        are invalid. All other exceptions lead to termination of the 
        optimization (see ASA IMMEDIATE_EXIT option [3]_).
    x0 : (N,) ndarray
        Initial parameters. This array is clobbered.
    xmin : (N,) ndarray
        Minimally allowed values for parameters
    xmin : (N,) ndarray
        Maximally allowed values for parameters
    full_output : bool, optional
        If True, return more information on completion,
        otherwise return only the final parameter values. See below.
    args : tuple, optional
        Additional positional arguments to pass to `func()`
    kwargs : dict, optional
        Additional keyword arguments to pass to `func()`
    parameter_type : (N,) ndarray, optional
        Array of type flags for each parameter:
            -1: REAL, -2: REAL_NO_REANNEAL
            1: INTEGER, 2: INTEGER_NO_REANNEAL

    Returns
    --------
    x : (N,) ndarray
        final parameter values
    f0 : float
        final cost value
    exit_code : int
        exit status code, look up in `asa_codes` to translate to human
        readable message.
    asa_opts : dict
        other parameters extracted from the ASA USER_DEFINES struct [3]_.

    Other Parameters
    ----------------
    All optional. See [3]_ for an explanation.

    rand_seed : long : 696969
    limit_acceptances : int : 1000
    limit_generated : int : 99999
    limit_invalid_generated_states : int : 1000
    accepted_to_generated_ratio : double : 1e-4
    cost_precision : double : 1e-18
    maximum_cost_repeat : int : 5
    number_cost_samples : int : 5
    temperature_ratio_scale : double : 1e-5
    cost_parameter_scale_ratio : double : 1.
    temperature_anneal_scale : double : 100. 
    include_integer_parameters : int : False
    user_initial_parameters : int : False
    sequential_parameters : int : -1
    initial_parameter_temperature : double : 1.
    acceptance_frequency_modulus : int : 100
    generated_frequency_modulus : int : 10000
    reanneal_cost : int : 1
    reanneal_parameters : int : 1
    delta_x : double : 1e-3

    Notes
    -----
    This algorithm permits an annealing schedule for "temperature" T
    decreasing exponentially in annealing-time k, T = T_0 exp(-c k^1/D).
    The introduction of re-annealing also permits adaptation to changing
    sensitivities in the multi-dimensional parameter-space. This
    annealing schedule is faster than fast Cauchy annealing, where T =
    T_0/k, and much faster than Boltzmann annealing, where T = T_0/ln k
    [2]_.

    References
    ----------
    For a description of the algorithm, documentation of options
    and tuning, see:

    .. [1] L. Ingber, Adaptive Simulated Annealing (ASA),
       Global optimization C-code, Caltech Alumni Association,
       Pasadena, CA, 1993, http://www.ingber.com/#ASA-CODE
    .. [2] L. Ingber, Very fast simulated re-annealing, 
       Mathematical Computer Modelling, vol 12, 8, pp 967-973, 1989,
       http://www.ingber.com/asa89_vfsr.pdf
    .. [3] L. Ingber, ASA Readme, http://www.ingber.com/ASA-README.pdf
    .. [4] A. Corana, M. Marchesi, C. Martini, S. Ridella,
       Minimizing multimodal functions of continuous variables with 
       the "simulated annealing" algorithm, ACM Trans. Mathl. Software,
       vol 13, 3, pp 262-279, 1987

    Examples
    --------
    The ASA test example [4]_:
    >>> s, t, c = .2, .05, .15
    >>> d = np.array([1., 1000., 10., 100.])
    >>> def cost(x):
    ...     k = np.rint(x/s)
    ...     r = np.fabs(k*s-x)
    ...     p = np.sign(k)
    ...     q = np.where(r<t, c*(p*p*k*s-p*t)**2, x**2)
    ...     return (d*q).sum()
    >>> x0 = np.array([999., -1007, 1001, -903])
    >>> xmax = 1e4*np.ones_like(x0)
    >>> x, f, err, opts = asa(cost, x0, -xmax, xmax, full_output=True)
    >>> print x, f, err # optimal parameters, cost, asa_code
    [ 0.01467097  0.02876664  0.02477989 -0.03838191] 0.0 2
    >>> print opts["n_accepted"], opts["n_generated"]
    335 3680
    """
    cdef USER_DEFINES opts
    cdef LONG_INT seed=rand_seed
    cdef int exit_code=0, cost_flag=1
    cdef int n=x0.shape[0]
    cdef ALLOC_INT param_num=n
    cdef np.ndarray curve=np.zeros([n, n], dtype=np.double)
    cdef np.ndarray tang=np.zeros([n], dtype=np.double)
    cdef double f0
    cdef void* data[4]

    data[0] = <void*>func
    data[1] = <void*>args
    data[2] = <void*>kwargs
    data[3] = NULL # exc_info here

    if parameter_type is None:
        parameter_type = -np.ones([n], dtype=np.int)

    opts.Limit_Acceptances = limit_acceptances
    opts.Limit_Generated = limit_generated
    opts.Limit_Invalid_Generated_States = limit_invalid_generated_states
    opts.Accepted_To_Generated_Ratio = accepted_to_generated_ratio

    opts.Cost_Precision = cost_precision
    opts.Maximum_Cost_Repeat = maximum_cost_repeat
    opts.Number_Cost_Samples = number_cost_samples
    opts.Temperature_Ratio_Scale = temperature_ratio_scale
    opts.Cost_Parameter_Scale_Ratio = cost_parameter_scale_ratio
    opts.Temperature_Anneal_Scale = temperature_anneal_scale

    opts.Include_Integer_Parameters = include_integer_parameters
    opts.User_Initial_Parameters = user_initial_parameters
    opts.Sequential_Parameters = sequential_parameters
    opts.Initial_Parameter_Temperature = initial_parameter_temperature

    opts.Acceptance_Frequency_Modulus = acceptance_frequency_modulus
    opts.Generated_Frequency_Modulus = generated_frequency_modulus
    opts.Reanneal_Cost = reanneal_cost
    opts.Reanneal_Parameters = reanneal_parameters

    opts.Delta_X = delta_x
    opts.User_Tangents = False
    opts.Curvature_0 = False
    #opts.Asa_Out_File = asa_out_file
    opts.Asa_Recursive_Level = 0

    opts.Asa_Data_Ptr = data
    opts.Asa_Data_Dim_Ptr = 1

    opts.Immediate_Exit = False

    resettable_randflt(&seed, 1)

    f0 = c_asa(cost_function, randflt, &seed,
            <double*>x0.data, <double*>xmin.data, <double*>xmax.data,
            <double*>tang.data, <double*>curve.data,
            &param_num, <int *>parameter_type.data,
            &cost_flag, &exit_code, &opts)

    asa_opts = dict(
            n_accepted=opts.N_Accepted,
            n_generated=opts.N_Generated,
            )

    if exit_code == IMMEDIATE_EXIT:
        if data[3] != NULL:
            exc_info = <object>data[3]
            Py_DECREF(exc_info)
            raise exc_info[0], exc_info[1], exc_info[2]
        raise SystemError("asa(): IMMEDIATE_EXIT without info")
    if full_output:
        return x0, f0, exit_code, asa_opts
    else:
        return x0


cdef double cost_function(double *x, double *xmin, double *xmax,
               double *tang, double *curve,
               ALLOC_INT *param_num, int *param_type,
               int *cost_flag, int *exit_code, USER_DEFINES * opts):
    cdef double cost
    cdef np.ndarray x_
    cdef object func, args, kwargs
    cdef np.npy_intp n = param_num[0]
    cdef void** data = <void**>opts.Asa_Data_Ptr

    if opts.Immediate_Exit:
        # c_asa() calls several times after we have set this
        return 0.

    func = <object>data[0]
    args = <object>data[1]
    kwargs = <object>data[2]
    x_ = np.PyArray_SimpleNewFromData(1, &n, np.NPY_DOUBLE, x)
    #xmin_ = np.PyArray_SimpleNewFromData(1, &n, np.NPY_DOUBLE, xmin)
    #xmax_ = np.PyArray_SimpleNewFromData(1, &n, np.NPY_DOUBLE, xmax)

    try:
        cost = func(x_, *args, **kwargs)
        cost_flag[0] = True
        return cost
    except CostParameterError:
        cost_flag[0] = False
        return 0.
    except:
        # catch all exceptions, save them and have asa() resurrect them
        exc_info = sys.exc_info()
        data[3] = <void*>exc_info
        Py_INCREF(exc_info)
        opts.Immediate_Exit = True
        return 0.
