/*
 * File: wls_alloc_mch.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 23-Nov-2019 12:19:19
 */

#ifndef WLS_ALLOC_MCH_H
#define WLS_ALLOC_MCH_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "wls_alloc_mch_types.h"

/* Function Declarations */
extern void wls_alloc_mch(const double v[3], double u[4], double p_limits,
  boolean_T v_limits);
extern void wls_alloc_mch_initialize(void);
extern void wls_alloc_mch_terminate(void);

#endif

/*
 * File trailer for wls_alloc_mch.h
 *
 * [EOF]
 */
