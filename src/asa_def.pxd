# -*- coding: utf8 -*-
#
#   pyasa - python bindings for Adaptive Simulated Annealing
#   Copyright (C) 2012 Robert Jordens <jordens@gmail.com>

cdef extern from "asa_usr_asa.h":
    ctypedef int BOOL
    ctypedef int ALLOC_INT
    ctypedef long LONG_INT
    ctypedef struct USER_DEFINES:
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
        char *Asa_Out_File
        #endif

        #if USER_ACCEPTANCE_TEST
        #void (*Acceptance_Test) (double cost,
        #                         double *parameter_minimum,
        #                         double *parameter_maximum,
        #                         ALLOC_INT * number_parameters,
        #                         const void *OPTIONS_TMP)
        #int User_Acceptance_Flag
        #int Cost_Acceptance_Flag
        #double Cost_Temp_Curr
        #double Cost_Temp_Init
        #double Cost_Temp_Scale
        #double Prob_Bias
        #LONG_INT *Random_Seed
        #endif

        #if ASA_SAMPLE
        #double Bias_Acceptance
        #double *Bias_Generated
        #double Average_Weights
        #double Limit_Weights
        #endif

        int Asa_Recursive_Level


cdef extern from "asa_rand.h":
    double randflt (LONG_INT * rand_seed)
    double resettable_randflt (LONG_INT * rand_seed, int reset)


cdef extern from "asa.h":
    double c_asa "asa" (double (*user_cost_function) (
            double *, double *, double *, double *, double *,
            ALLOC_INT *, int *, int *, int *, USER_DEFINES *),
         double (*user_random_generator) (LONG_INT *), LONG_INT * rand_seed,
         double *parameter_initial_final, double *parameter_minimum,
         double *parameter_maximum, double *tangents, double *curvature,
         ALLOC_INT * number_parameters, int *parameter_type,
         int *valid_state_generated_flag, int *exit_status,
         USER_DEFINES * OPTIONS)
    
    cdef enum ASE_EXIT_CODES:
        NORMAL_EXIT = 0
        P_TEMP_TOO_SMALL = 1
        C_TEMP_TOO_SMALL = 2
        COST_REPEATING = 3
        TOO_MANY_INVALID_STATES = 4
        IMMEDIATE_EXIT = 5
        INVALID_USER_INPUT = 7
        INVALID_COST_FUNCTION = 8
        INVALID_COST_FUNCTION_DERIV = 9
        CALLOC_FAILED = -1
