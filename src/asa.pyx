# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Lester Ingber's Adaptive Stimulated
#   Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@phys.ethz.ch>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
cimport numpy as np

from cpython cimport PyObject, Py_INCREF

np.import_array()

cdef extern from "asa_usr_asa.h":
    ctypedef int BOOL
    ctypedef int ALLOC_INT
    ctypedef long LONG_INT
    cdef struct _USER_DEFINES:
        LONG_INT Limit_Acceptances
        LONG_INT Limit_Generated
        int Limit_Invalid_Generated_States
        double Accepted_To_Generated_Ratio

        double Cost_Precision
        int Maximum_Cost_Repeat
        int Number_Cost_Samples
        double Temperature_Ratio_Scale
        double Cost_Parameter_Scale_Ratio
        double Temperature_Anneal_Scale

        int Include_Integer_Parameters
        int User_Initial_Parameters
        ALLOC_INT Sequential_Parameters
        double Initial_Parameter_Temperature

        int Acceptance_Frequency_Modulus
        int Generated_Frequency_Modulus
        int Reanneal_Cost
        int Reanneal_Parameters

        double Delta_X
        int User_Tangents
        int Curvature_0

        #if QUENCH_PARAMETERS
        #double *User_Quench_Param_Scale
        #endif
        #if QUENCH_COST
        #double *User_Quench_Cost_Scale
        #endif

        LONG_INT N_Accepted
        LONG_INT N_Generated
        int Locate_Cost
        int Immediate_Exit

        double *Best_Cost
        double *Best_Parameters
        double *Last_Cost
        double *Last_Parameters

        ALLOC_INT Asa_Data_Dim_Ptr
        void * Asa_Data_Ptr
        #if USER_ASA_OUT
        #char *Asa_Out_File
        #endif
        int Asa_Recursive_Level

    ctypedef _USER_DEFINES USER_DEFINES


cdef extern from "asa_rand.h":
    double randflt (LONG_INT * rand_seed)


cdef extern from "asa.h":
    double c_asa "asa" (double (*user_cost_function)
         (double *, double *, double *, double *, double *, ALLOC_INT *,
          int *, int *, int *, USER_DEFINES *),
         double (*user_random_generator) (LONG_INT *), LONG_INT * rand_seed,
         double *parameter_initial_final, double *parameter_minimum,
         double *parameter_maximum, double *tangents, double *curvature,
         ALLOC_INT * number_parameters, int *parameter_type,
         int *valid_state_generated_flag, int *exit_status,
         USER_DEFINES * OPTIONS)


def asa(func, x0, **kw):
    cdef USER_DEFINES OPTIONS
    cdef LONG_INT rand_seed = kw.get("seed", 0)
    cdef int exit_status, valid_state_gen_flag, param_type
    cdef ALLOC_INT num_param = len(x0)
    cdef double param_init_fin, param_min, param_max, tang, curve
    OPTIONS.Asa_Data_Ptr = <void*>func
    return c_asa(cost_function, randflt, &rand_seed,
            &param_init_fin, &param_min, &param_max,
            &tang, &curve,
            &num_param, &param_type,
            &valid_state_gen_flag, &exit_status, &OPTIONS)


cdef double cost_function (double *x,
               double *parameter_lower_bound,
               double *parameter_upper_bound,
               double *cost_tangents,
               double *cost_curvature,
               ALLOC_INT * parameter_dimension,
               int *parameter_int_real,
               int *cost_flag, int *exit_code,
               USER_DEFINES * USER_OPTIONS):
    return (<object>USER_OPTIONS.Asa_Data_Ptr)()
