# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Adaptive Simulated Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@gmail.com>

import sys

import numpy as np
cimport numpy as np

from cpython cimport Py_INCREF, Py_DECREF

from asa_def cimport *

np.import_array()

asa_errors = {

    NORMAL_EXIT: """Given the criteria set, the
search has run its normal course, e.g., when Limit_Acceptances or
Limit_Generated is reached.""", 

    P_TEMP_TOO_SMALL: """A parameter temperature was too small using the
set criteria.  Often this is an acceptable status code.""",

    C_TEMP_TOO_SMALL: """The cost temperature was too small using the set
criteria.  Often this is an acceptable status code.""",

    COST_REPEATING: """The cost function value repeated a number of times
using the set criteria. Often this is an acceptable status code.""",

    TOO_MANY_INVALID_STATES: """Too many repetitive generated states were
invalid using the set criteria. This is helpful when using *cost_flag, as
discussed above, to include constraints.""",

    IMMEDIATE_EXIT: """The user has set OPTIONS−>Immediate_Exit to TRUE, or
deleted file asa_exit_anytime when ASA_EXIT_ANYTIME is TRUE.""",

    INVALID_USER_INPUT: """The user has introduced invalid input. When
entering asa (), a function asa_test_asa_options () checks out many
user−defined parameters and OPTIONS, and prints out invalid OPTIONS when
ASA_PRINT is set to TRUE.""",

    INVALID_COST_FUNCTION: """The user has returned a value of the cost
function to asa () which is not a valid number, e.g., not between
-MAX_DOUBLE and MAX_DOUBLE.  Or, the user has returned a value of a
parameter no longer within its proper range (excluding cases where the
user has set the lower bound equal to the upper bound to remove a parameter
from consideration).""",

    INVALID_COST_FUNCTION_DERIV: """While calculating numerical cost
derivatives, a value of the cost function was returned which is not a
valid number, e.g., not between -MAX_DOUBLE and MAX_DOUBLE. Or, while
calculating numerical cost derivatives, a value of a parameter no longer
within its proper range (excluding cases where the user has set the
lower bound equal to the upper bound to remove a parameter from
consideration) was set.  In such cases, review the bounds of parameters
and the OPTIONS used to determine how derivatives are calculated and
used.  When calculating derivatives, an invalid cost function, i.e.,
returning *cost_flag = FALSE from the cost function, will exit asa()
with *exit_code INVALID_COST_FUNCTION_DERIV. Also, when calculating
derivatives, no extra test is performed to check that parameters are
within their lower and upper bounds (since meshes for derivatives may
exceed these bounds). If this is not desired, then within the cost
function, within a test for USER_OPTIONS->Locate_Cost == 3, a decision
may be made whether to return this *exit_code.""",

    CALLOC_FAILED: """Calloc memory allocation has failed in asa.c. This error
will call Exit_ASA(), the location will be printed to stdout, and asa () will
return the double -1 to the calling program. In asa_usr.c, if asa () returns
this *exit_code a warning will be printed both to stdout and to USER_OUT.
Note that if a calloc memory allocation fails in asa_usr.c, this error will
call Exit_USER() to print the location to stdout and then return -2.""",

}


class CostParameterError(Exception):
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
    http://www.ingber.com/#ASA
        
    Parameters:
    -----------

        object func not None,
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

    Returns:
    --------

        x0, f0, exit_code, asa_errors[exit_code], curve, asa_opts
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
    opts.Asa_Data_Dim_Ptr = 4

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
        else:
            raise SystemError("asa(): IMMEDIATE_EXIT without info")
    if full_output:
        return x0, f0, exit_code, asa_errors[exit_code], curve, asa_opts
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
