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
    cdef struct USER_DEFINES:
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
        #if USER_INITIAL_COST_TEMP
        #double *User_Cost_Temperature
        ##endif

        int Include_Integer_Parameters
        int User_Initial_Parameters
        ALLOC_INT Sequential_Parameters
        double Initial_Parameter_Temperature
        #if RATIO_TEMPERATURE_SCALES
        #double *User_Temperature_Ratio
        #endif
        #if USER_INITIAL_PARAMETERS_TEMPS
        #double *User_Parameter_Temperature
        #endif

        int Acceptance_Frequency_Modulus
        int Generated_Frequency_Modulus
        int Reanneal_Cost
        int Reanneal_Parameters

        double Delta_X
        #if DELTA_PARAMETERS
        #double *User_Delta_Parameter
        #endif
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

        #if OPTIONAL_DATA_DBL
        #ALLOC_INT Asa_Data_Dim_Dbl
        #double *Asa_Data_Dbl
        #endif
        #if OPTIONAL_DATA_INT
        #ALLOC_INT Asa_Data_Dim_Int
        #LONG_INT *Asa_Data_Int
        #endif
        #if OPTIONAL_DATA_PTR
        #ALLOC_INT Asa_Data_Dim_Ptr
        #OPTIONAL_PTR_TYPE *Asa_Data_Ptr
        #endif
        #if USER_ASA_OUT
        #char *Asa_Out_File
        #endif
        #if USER_COST_SCHEDULE
        #if HAVE_ANSI
        #double (*Cost_Schedule) (double current_cost_temperature,
        #                         const void *OPTIONS_TMP)
        #else                           /* HAVE_ANSI */
        #double (*Cost_Schedule) ()
        #endif                          /* HAVE_ANSI */
        #endif
        #if USER_ACCEPT_ASYMP_EXP
        #double Asymp_Exp_Param
        #endif
        #if USER_ACCEPTANCE_TEST
        #if HAVE_ANSI
        #void (*Acceptance_Test) (double cost,
        #                         double *parameter_minimum,
        #                         double *parameter_maximum,
        #                         ALLOC_INT * number_parameters,
        #                         const void *OPTIONS_TMP)
        #else                           /* HAVE_ANSI */
        #void (*Acceptance_Test) ()
        ##endif                          /* HAVE_ANSI */
        #int User_Acceptance_Flag
        #int Cost_Acceptance_Flag
        #double Cost_Temp_Curr
        #double Cost_Temp_Init
        #double Cost_Temp_Scale
        #double Prob_Bias
        #LONG_INT *Random_Seed
        #endif
        #if USER_GENERATING_FUNCTION
        #if HAVE_ANSI
        #double (*Generating_Distrib) (LONG_INT * seed,
        #                              ALLOC_INT * parameter_dimension,
        #                              ALLOC_INT index_v,
        #                              double temperature_v,
        #                              double init_param_temp_v,
        #                              double temp_scale_params_v,
        #                              double parameter_v,
        #                              double parameter_range_v,
        #                              double *last_saved_parameter,
        #                              const void *OPTIONS_TMP)
        #else                           /* HAVE_ANSI */
        #double (*Generating_Distrib) ()
        #endif                          /* HAVE_ANSI */
        #endif
        #if USER_REANNEAL_COST
        #if HAVE_ANSI
        #int (*Reanneal_Cost_Function) (double *cost_best,
        #                               double *cost_last,
        #                               double *initial_cost_temperature,
        #                               double *current_cost_temperature,
        #                               const void *OPTIONS_TMP)
        #else                           /* HAVE_ANSI */
        #int (*Reanneal_Cost_Function) ()
        #endif                          /* HAVE_ANSI */
        #endif
        #if USER_REANNEAL_PARAMETERS
        #if HAVE_ANSI
        #double (*Reanneal_Params_Function) (double current_temp,
        #                                    double tangent,
        #                                    double max_tangent,
        #                                    const void *OPTIONS_TMP)
        #else                           /* HAVE_ANSI */
        #double (*Reanneal_Params_Function) ()
        #endif                          /* HAVE_ANSI */
        #endif
        #if ASA_SAMPLE
        #double Bias_Acceptance
        #double *Bias_Generated
        #double Average_Weights
        #double Limit_Weights
        #endif
        #if ASA_QUEUE
        #ALLOC_INT Queue_Size
        #double *Queue_Resolution
        #endif
        #if ASA_RESOLUTION
        #double *Coarse_Resolution
        #endif
        #if FITLOC
        #int Fit_Local
        #int Iter_Max
        #double Penalty
        #endif
        #if MULTI_MIN
        #int Multi_Number
        #double *Multi_Cost
        #double **Multi_Params
        #double *Multi_Grid
        #int Multi_Specify
        #endif
        #if ASA_PARALLEL
        #int Gener_Mov_Avr
        #LONG_INT Gener_Block
        #LONG_INT Gener_Block_Max
        #endif
        #if ASA_SAVE
        #ALLOC_INT Random_Array_Dim
        #double *Random_Array
        #endif
        int Asa_Recursive_Level
        #if ASA_FUZZY
        #int NoOfSamples
        #double ThresholdDeviation
        #double Threshold1
        #double Performance_Target
        #double Factor_a
        #endif

cdef extern from "asa.h":
    double asa (double (*user_cost_function)
         (double *, double *, double *, double *, double *, ALLOC_INT *,
          int *, int *, int *, USER_DEFINES *),
         double (*user_random_generator) (LONG_INT *), LONG_INT * rand_seed,
         double *parameter_initial_final, double *parameter_minimum,
         double *parameter_maximum, double *tangents, double *curvature,
         ALLOC_INT * number_parameters, int *parameter_type,
         int *valid_state_generated_flag, int *exit_status,
         USER_DEFINES * OPTIONS)


