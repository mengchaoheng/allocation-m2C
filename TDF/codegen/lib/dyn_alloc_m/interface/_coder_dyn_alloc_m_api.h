/*
 * File: _coder_dyn_alloc_m_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 18-Aug-2020 19:12:23
 */

#ifndef _CODER_DYN_ALLOC_M_API_H
#define _CODER_DYN_ALLOC_M_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_dyn_alloc_m_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void dyn_alloc_m(real_T B[27], real_T v[3], real_T u[9], real_T p_limits
  [2], real_T v_limits[2], real_T T, real_T Wv[9], real_T W1[81], real_T W2[81],
  real_T S[27], real_T imax, real_T gam, boolean_T only_plim);
extern void dyn_alloc_m_api(const mxArray * const prhs[13], int32_T nlhs, const
  mxArray *plhs[1]);
extern void dyn_alloc_m_atexit(void);
extern void dyn_alloc_m_initialize(void);
extern void dyn_alloc_m_terminate(void);
extern void dyn_alloc_m_xil_terminate(void);

#endif

/*
 * File trailer for _coder_dyn_alloc_m_api.h
 *
 * [EOF]
 */
