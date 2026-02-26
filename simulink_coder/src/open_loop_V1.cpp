/*
 * open_loop_V1.cpp
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "open_loop_V1".
 *
 * Model version              : 12.20
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Thu Feb 26 18:13:08 2026
 *
 * Target selection: ert.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "open_loop_V1.h"
#include "open_loop_V1_types.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include <emmintrin.h>
#include "rmw/qos_profiles.h"
#include <stddef.h>

extern "C"
{

#include "rt_nonfinite.h"

}

#include "rt_defines.h"
#include "open_loop_V1_private.h"

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
void open_loop_V1::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = static_cast<ODE3_IntgData *>(rtsiGetSolverData(si));
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 17;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                static_cast<uint_T>(nXc)*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  open_loop_V1_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  open_loop_V1_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  this->step();
  open_loop_V1_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

void open_loop_V1::open_lo_ServiceCaller_setupImpl(const
  ros_slros2_internal_block_Ser_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[25] = "/gazebo/set_entity_state";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S2>/ServiceCaller' */
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 25; i++) {
    /* Start for MATLABSystem: '<S2>/ServiceCaller' */
    open_loop_V1_B.b_zeroDelimTopic[i] = b_zeroDelimTopic[i];
  }

  ServCall_open_loop_V1_326.createServiceCaller
    (&open_loop_V1_B.b_zeroDelimTopic[0], qos_profile);
}

real_T open_loop_V1::open_loop_V1_rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      open_loop_V1_B.i_g = 1;
    } else {
      open_loop_V1_B.i_g = -1;
    }

    if (u1 > 0.0) {
      open_loop_V1_B.i1 = 1;
    } else {
      open_loop_V1_B.i1 = -1;
    }

    y = atan2(static_cast<real_T>(open_loop_V1_B.i_g), static_cast<real_T>
              (open_loop_V1_B.i1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T open_loop_V1::open_loop_V1_rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    open_loop_V1_B.d = fabs(u0);
    open_loop_V1_B.d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (open_loop_V1_B.d == 1.0) {
        y = 1.0;
      } else if (open_loop_V1_B.d > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (open_loop_V1_B.d1 == 0.0) {
      y = 1.0;
    } else if (open_loop_V1_B.d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Model step function */
void open_loop_V1::step()
{
  static const real_T a[9] = { 0.0256247963, 0.0, -0.0020093863, 0.0,
    0.0116945386, 0.0, -0.0020093863, 0.0, 0.0088954752 };

  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  static const uint8_T b[11] = { 101U, 107U, 114U, 97U, 110U, 111U, 112U, 108U,
    97U, 110U, 111U };

  static const uint8_T b_0[5] = { 119U, 111U, 114U, 108U, 100U };

  static const real_T a_0[9] = { 0.0256247963, 0.0, -0.0020093863, 0.0,
    0.0116945386, 0.0, -0.0020093863, 0.0, 0.0088954752 };

  static const int8_T c[3] = { 1, 0, 0 };

  if (rtmIsMajorTimeStep((&open_loop_V1_M))) {
    /* set solver stop time */
    if (!((&open_loop_V1_M)->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&(&open_loop_V1_M)->solverInfo, (((&open_loop_V1_M
        )->Timing.clockTickH0 + 1) * (&open_loop_V1_M)->Timing.stepSize0 *
        4294967296.0));
    } else {
      rtsiSetSolverStopTime(&(&open_loop_V1_M)->solverInfo, (((&open_loop_V1_M
        )->Timing.clockTick0 + 1) * (&open_loop_V1_M)->Timing.stepSize0 +
        (&open_loop_V1_M)->Timing.clockTickH0 * (&open_loop_V1_M)
        ->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&open_loop_V1_M))) {
    (&open_loop_V1_M)->Timing.t[0] = rtsiGetT(&(&open_loop_V1_M)->solverInfo);
  }

  open_loop_V1_B.b = rtmIsMajorTimeStep((&open_loop_V1_M));
  if (open_loop_V1_B.b) {
    /* MATLAB Function: '<Root>/MATLAB Function' */
    memset(&open_loop_V1_B.stringOut_l[0], 0, sizeof(uint8_T) << 7U);
    for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 11; open_loop_V1_B.i++) {
      open_loop_V1_B.stringOut_l[open_loop_V1_B.i] = b[open_loop_V1_B.i];
    }

    open_loop_V1_B.lengthOut_e = 11U;

    /* End of MATLAB Function: '<Root>/MATLAB Function' */

    /* MATLAB Function: '<Root>/MATLAB Function1' */
    memset(&open_loop_V1_B.stringOut[0], 0, sizeof(uint8_T) << 7U);
    for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 5; open_loop_V1_B.i++) {
      open_loop_V1_B.stringOut[open_loop_V1_B.i] = b_0[open_loop_V1_B.i];
    }

    open_loop_V1_B.lengthOut = 5U;

    /* End of MATLAB Function: '<Root>/MATLAB Function1' */
  }

  /* Integrator: '<S7>/Integrator' */
  memcpy(&open_loop_V1_B.x[0], &open_loop_V1_X.Integrator_CSTATE[0], 12U *
         sizeof(real_T));

  /* Gain: '<S97>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Gain: '<S87>/Derivative Gain'
   *  Integrator: '<S89>/Filter'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<S89>/SumD'
   */
  open_loop_V1_B.FilterCoefficient = ((1.0 - open_loop_V1_B.x[11]) * 0.0001 -
    open_loop_V1_X.Filter_CSTATE) * 100.0;

  /* Sum: '<S103>/Sum' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Gain: '<S99>/Proportional Gain'
   *  Integrator: '<S94>/Integrator'
   *  Sum: '<Root>/Sum2'
   */
  open_loop_V1_B.Saturation = ((1.0 - open_loop_V1_B.x[11]) * 1.0E-5 +
    open_loop_V1_X.Integrator_CSTATE_b) + open_loop_V1_B.FilterCoefficient;

  /* Saturate: '<S101>/Saturation' */
  if (open_loop_V1_B.Saturation > 0.000125) {
    /* Sum: '<S103>/Sum' incorporates:
     *  Saturate: '<S101>/Saturation'
     */
    open_loop_V1_B.Saturation = 0.000125;
  } else if (open_loop_V1_B.Saturation < 0.0) {
    /* Sum: '<S103>/Sum' incorporates:
     *  Saturate: '<S101>/Saturation'
     */
    open_loop_V1_B.Saturation = 0.0;
  }

  /* End of Saturate: '<S101>/Saturation' */

  /* Gain: '<S45>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S35>/Derivative Gain'
   *  Integrator: '<S37>/Filter'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S37>/SumD'
   */
  open_loop_V1_B.FilterCoefficient_m = ((-0.034906585039886591 -
    open_loop_V1_B.x[7]) * 0.5 - open_loop_V1_X.Filter_CSTATE_m) * 100.0;

  /* Sum: '<S51>/Sum' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S47>/Proportional Gain'
   *  Integrator: '<S42>/Integrator'
   *  Sum: '<Root>/Sum1'
   */
  open_loop_V1_B.Saturation_f = ((-0.034906585039886591 - open_loop_V1_B.x[7]) *
    2.5 + open_loop_V1_X.Integrator_CSTATE_p) +
    open_loop_V1_B.FilterCoefficient_m;

  /* Saturate: '<S49>/Saturation' */
  if (open_loop_V1_B.Saturation_f > 0.3490658503988659) {
    /* Sum: '<S51>/Sum' incorporates:
     *  Saturate: '<S49>/Saturation'
     */
    open_loop_V1_B.Saturation_f = 0.3490658503988659;
  } else if (open_loop_V1_B.Saturation_f < -0.3490658503988659) {
    /* Sum: '<S51>/Sum' incorporates:
     *  Saturate: '<S49>/Saturation'
     */
    open_loop_V1_B.Saturation_f = -0.3490658503988659;
  }

  /* End of Saturate: '<S49>/Saturation' */

  /* SignalConversion generated from: '<S112>/ SFunction ' incorporates:
   *  Constant: '<Root>/nominal_control'
   *  MATLAB Function: '<S7>/MATLAB Function'
   */
  open_loop_V1_B.control_vector[0] = 0.0;
  open_loop_V1_B.control_vector[1] = open_loop_V1_B.Saturation_f;
  open_loop_V1_B.control_vector[2] = 0.0;
  open_loop_V1_B.control_vector[3] = open_loop_V1_B.Saturation;
  open_loop_V1_B.control_vector[4] = open_loop_V1_B.Saturation;

  /* MATLAB Function: '<S7>/MATLAB Function' */
  open_loop_V1_B.u2 = open_loop_V1_B.control_vector[1];
  open_loop_V1_B.u4 = open_loop_V1_B.control_vector[3];
  open_loop_V1_B.u5 = open_loop_V1_B.control_vector[4];
  if (open_loop_V1_B.control_vector[1] > 0.3490658503988659) {
    open_loop_V1_B.u2 = 0.3490658503988659;
  } else if (open_loop_V1_B.control_vector[1] < -0.3490658503988659) {
    open_loop_V1_B.u2 = -0.3490658503988659;
  }

  if (open_loop_V1_B.control_vector[3] > 0.000125) {
    open_loop_V1_B.u4 = 0.000125;
  } else if (open_loop_V1_B.control_vector[3] < 0.0) {
    open_loop_V1_B.u4 = 0.0;
  }

  if (open_loop_V1_B.control_vector[4] > 0.000125) {
    open_loop_V1_B.u5 = 0.000125;
  } else if (open_loop_V1_B.control_vector[4] < 0.0) {
    open_loop_V1_B.u5 = 0.0;
  }

  open_loop_V1_B.c_phi = cos(open_loop_V1_B.x[6]);
  open_loop_V1_B.s_phi = sin(open_loop_V1_B.x[6]);
  open_loop_V1_B.c_theta = cos(open_loop_V1_B.x[7]);
  open_loop_V1_B.s_theta = sin(open_loop_V1_B.x[7]);
  open_loop_V1_B.c_psi = cos(open_loop_V1_B.x[8]);
  open_loop_V1_B.s_psi = sin(open_loop_V1_B.x[8]);
  open_loop_V1_B.c_theta_tmp_f = open_loop_V1_B.c_theta * open_loop_V1_B.c_psi;
  open_loop_V1_B.c_theta_g[0] = open_loop_V1_B.c_theta_tmp_f;
  open_loop_V1_B.c_theta_tmp_cv = open_loop_V1_B.c_theta * open_loop_V1_B.s_psi;
  open_loop_V1_B.c_theta_g[3] = open_loop_V1_B.c_theta_tmp_cv;
  open_loop_V1_B.c_theta_g[6] = -open_loop_V1_B.s_theta;
  open_loop_V1_B.c_theta_tmp = open_loop_V1_B.s_phi * open_loop_V1_B.s_theta;
  open_loop_V1_B.c_theta_tmp_p = open_loop_V1_B.c_theta_tmp *
    open_loop_V1_B.c_psi - open_loop_V1_B.c_phi * open_loop_V1_B.s_psi;
  open_loop_V1_B.c_theta_g[1] = open_loop_V1_B.c_theta_tmp_p;
  open_loop_V1_B.c_theta_tmp = open_loop_V1_B.c_theta_tmp * open_loop_V1_B.s_psi
    + open_loop_V1_B.c_phi * open_loop_V1_B.c_psi;
  open_loop_V1_B.c_theta_g[4] = open_loop_V1_B.c_theta_tmp;
  open_loop_V1_B.c_theta_tmp_b = open_loop_V1_B.s_phi * open_loop_V1_B.c_theta;
  open_loop_V1_B.c_theta_g[7] = open_loop_V1_B.c_theta_tmp_b;
  open_loop_V1_B.c_theta_tmp_c = open_loop_V1_B.c_phi * open_loop_V1_B.s_theta;
  open_loop_V1_B.c_theta_tmp_cx = open_loop_V1_B.c_theta_tmp_c *
    open_loop_V1_B.c_psi + open_loop_V1_B.s_phi * open_loop_V1_B.s_psi;
  open_loop_V1_B.c_theta_g[2] = open_loop_V1_B.c_theta_tmp_cx;
  open_loop_V1_B.c_theta_tmp_c = open_loop_V1_B.c_theta_tmp_c *
    open_loop_V1_B.s_psi - open_loop_V1_B.s_phi * open_loop_V1_B.c_psi;
  open_loop_V1_B.c_theta_g[5] = open_loop_V1_B.c_theta_tmp_c;
  open_loop_V1_B.c_theta_tmp_k = open_loop_V1_B.c_phi * open_loop_V1_B.c_theta;
  open_loop_V1_B.c_theta_g[8] = open_loop_V1_B.c_theta_tmp_k;
  for (open_loop_V1_B.i = 0; open_loop_V1_B.i <= 0; open_loop_V1_B.i += 2) {
    tmp = _mm_loadu_pd(&open_loop_V1_B.c_theta_g[open_loop_V1_B.i + 3]);
    tmp_2 = _mm_set1_pd(0.0);
    tmp_0 = _mm_loadu_pd(&open_loop_V1_B.c_theta_g[open_loop_V1_B.i]);
    tmp_1 = _mm_loadu_pd(&open_loop_V1_B.c_theta_g[open_loop_V1_B.i + 6]);
    _mm_storeu_pd(&open_loop_V1_B.Vb_w[open_loop_V1_B.i], _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp, tmp_2), _mm_mul_pd(tmp_0, tmp_2)), _mm_mul_pd(tmp_1,
      tmp_2)));
  }

  for (open_loop_V1_B.i = 2; open_loop_V1_B.i < 3; open_loop_V1_B.i++) {
    open_loop_V1_B.Vb_w[open_loop_V1_B.i] =
      (open_loop_V1_B.c_theta_g[open_loop_V1_B.i + 3] * 0.0 +
       open_loop_V1_B.c_theta_g[open_loop_V1_B.i] * 0.0) +
      open_loop_V1_B.c_theta_g[open_loop_V1_B.i + 6] * 0.0;
  }

  tmp = _mm_sub_pd(_mm_loadu_pd(&open_loop_V1_B.x[0]), _mm_loadu_pd
                   (&open_loop_V1_B.Vb_w[0]));
  _mm_storeu_pd(&open_loop_V1_B.dv[0], tmp);

  /* MATLAB Function: '<S7>/MATLAB Function' */
  open_loop_V1_B.Va_b_idx_2 = open_loop_V1_B.x[2] - open_loop_V1_B.Vb_w[2];
  open_loop_V1_B.Va = sqrt((open_loop_V1_B.dv[0] * open_loop_V1_B.dv[0] +
    open_loop_V1_B.dv[1] * open_loop_V1_B.dv[1]) + open_loop_V1_B.Va_b_idx_2 *
    open_loop_V1_B.Va_b_idx_2);
  open_loop_V1_B.c_psi = open_loop_V1_rt_atan2d_snf(open_loop_V1_B.Va_b_idx_2,
    open_loop_V1_B.dv[0]);
  open_loop_V1_B.s_psi = asin(open_loop_V1_B.dv[1] / open_loop_V1_B.Va);
  open_loop_V1_B.Va_b_idx_2 = open_loop_V1_B.Va * open_loop_V1_B.Va * 0.6125;
  open_loop_V1_B.wbe_b[0] = open_loop_V1_B.x[3];
  open_loop_V1_B.wbe_b[1] = open_loop_V1_B.x[4];
  open_loop_V1_B.wbe_b[2] = open_loop_V1_B.x[5];
  open_loop_V1_B.CL_w_OGE = ((open_loop_V1_B.c_psi - -0.065449846949787352) +
    0.082903139469730644) * 4.9604094530365153;
  open_loop_V1_B.CL_h_OGE = ((open_loop_V1_B.c_psi - -0.074176493209759012) +
    0.043633231299858237) * 4.8387748917360032;
  open_loop_V1_B.CD_iw_IGE = fabs((-open_loop_V1_B.x[11] - 0.363) / 5.02);
  open_loop_V1_B.CL_w_IGE = (288.0 * open_loop_V1_rt_powd_snf
    (open_loop_V1_B.CD_iw_IGE, 0.787) * exp(-9.14 * open_loop_V1_rt_powd_snf
    (open_loop_V1_B.CD_iw_IGE, 0.327)) * 0.97986308862072491 /
    5.9129476540958859 + 1.0) * open_loop_V1_B.CL_w_OGE;
  open_loop_V1_B.CD_ih_IGE = fabs((-open_loop_V1_B.x[11] + 0.72) / 2.74);
  open_loop_V1_B.u2 = (288.0 * open_loop_V1_rt_powd_snf(open_loop_V1_B.CD_ih_IGE,
    0.787) * exp(-9.14 * open_loop_V1_rt_powd_snf(open_loop_V1_B.CD_ih_IGE,
    0.327)) * 0.95628590200128227 / 5.35300902982722 + 1.0) *
    ((open_loop_V1_B.u2 * open_loop_V1_B.u2 * -0.00141 + 0.0307 *
      open_loop_V1_B.u2) + open_loop_V1_B.CL_h_OGE);
  open_loop_V1_B.CD_iw_IGE = open_loop_V1_B.CL_w_OGE * open_loop_V1_B.CL_w_OGE /
    21.205750411731103 * (1.0 - exp(-10.1 * open_loop_V1_rt_powd_snf
    (open_loop_V1_B.CD_iw_IGE, 0.686)));
  open_loop_V1_B.CD_ih_IGE = open_loop_V1_B.CL_h_OGE * open_loop_V1_B.CL_h_OGE /
    18.943803701146454 * (1.0 - exp(-10.1 * open_loop_V1_rt_powd_snf
    (open_loop_V1_B.CD_ih_IGE, 0.686)));
  open_loop_V1_B.FA_b_tmp = sin(open_loop_V1_B.c_psi);
  open_loop_V1_B.FA_b_tmp_m = cos(open_loop_V1_B.c_psi);
  open_loop_V1_B.c_theta_g[0] = open_loop_V1_B.FA_b_tmp_m;
  open_loop_V1_B.c_theta_g[3] = 0.0;
  open_loop_V1_B.c_theta_g[6] = -open_loop_V1_B.FA_b_tmp;
  open_loop_V1_B.c_theta_g[2] = open_loop_V1_B.FA_b_tmp;
  open_loop_V1_B.c_theta_g[5] = 0.0;
  open_loop_V1_B.c_theta_g[8] = open_loop_V1_B.FA_b_tmp_m;
  open_loop_V1_B.CD_iw_IGE_m[0] = -(((open_loop_V1_B.CD_iw_IGE + 0.0306) +
    (open_loop_V1_B.CD_ih_IGE + 0.0306)) * open_loop_V1_B.Va_b_idx_2 * 3.334);
  open_loop_V1_B.CD_iw_IGE_m[1] = 0.0 * open_loop_V1_B.Va_b_idx_2 * 3.334;
  open_loop_V1_B.CD_iw_IGE_m[2] = -((open_loop_V1_B.CL_w_IGE * 3.334 +
    open_loop_V1_B.u2 * 1.128) * open_loop_V1_B.Va_b_idx_2);
  open_loop_V1_B.c_theta_g[1] = 0.0;
  open_loop_V1_B.FA_b_idx_0 = 0.0;
  open_loop_V1_B.c_theta_g[4] = 1.0;
  open_loop_V1_B.FA_b_idx_1 = 0.0;
  open_loop_V1_B.c_theta_g[7] = 0.0;
  open_loop_V1_B.FA_b_idx_2 = 0.0;
  for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 3; open_loop_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&open_loop_V1_B.c_theta_g[3 *
      open_loop_V1_B.i]), _mm_set1_pd
      (open_loop_V1_B.CD_iw_IGE_m[open_loop_V1_B.i])), _mm_set_pd
                     (open_loop_V1_B.FA_b_idx_1, open_loop_V1_B.FA_b_idx_0));
    _mm_storeu_pd(&open_loop_V1_B.dv[0], tmp);
    open_loop_V1_B.FA_b_idx_0 = open_loop_V1_B.dv[0];
    open_loop_V1_B.FA_b_idx_1 = open_loop_V1_B.dv[1];
    open_loop_V1_B.FA_b_idx_2 += open_loop_V1_B.c_theta_g[3 * open_loop_V1_B.i +
      2] * open_loop_V1_B.CD_iw_IGE_m[open_loop_V1_B.i];
  }

  open_loop_V1_B.Cl = -0.0005 * open_loop_V1_B.s_psi * 180.0 /
    3.1415926535897931;
  open_loop_V1_B.Cm = -0.02 * open_loop_V1_B.c_psi * 180.0 / 3.1415926535897931;
  open_loop_V1_B.Cn = -0.002 * open_loop_V1_B.s_psi * 180.0 / 3.1415926535897931;
  open_loop_V1_B.u4 = open_loop_V1_B.u4 / 0.000125 * (37.42 - open_loop_V1_B.Va)
    + open_loop_V1_B.Va;
  open_loop_V1_B.u5 = open_loop_V1_B.u5 / 0.000125 * (37.42 - open_loop_V1_B.Va)
    + open_loop_V1_B.Va;
  open_loop_V1_B.Vb_w[0] = 0.11056515539994617 * open_loop_V1_B.u4 *
    (open_loop_V1_B.u4 - open_loop_V1_B.Va) * 0.99619469809174555;
  open_loop_V1_B.u5 = 0.11056515539994617 * open_loop_V1_B.u5 *
    (open_loop_V1_B.u5 - open_loop_V1_B.Va) * 0.99619469809174555;
  open_loop_V1_B.u4 = 0.0;
  open_loop_V1_B.Va = -9.81 * open_loop_V1_B.s_theta * 112.0;
  open_loop_V1_B.c_theta *= 9.81;
  open_loop_V1_B.s_phi = open_loop_V1_B.c_theta * open_loop_V1_B.s_phi * 112.0;
  open_loop_V1_B.c_phi = open_loop_V1_B.c_theta * open_loop_V1_B.c_phi * 112.0;
  open_loop_V1_B.Mcg_b_idx_2 = 5.02 * open_loop_V1_B.Va_b_idx_2 * 3.334;
  open_loop_V1_B.c_theta = open_loop_V1_B.Mcg_b_idx_2 * open_loop_V1_B.Cl;
  open_loop_V1_B.Mcg_b_idx_1 = 0.646 * open_loop_V1_B.Va_b_idx_2 * 3.334 *
    open_loop_V1_B.Cm + ((-0.285 * open_loop_V1_B.Vb_w[0] - 0.04523383048603459)
    + (-0.285 * open_loop_V1_B.u5 - 0.04523383048603459));
  open_loop_V1_B.Mcg_b_idx_2 = ((0.0 - 0.6 * open_loop_V1_B.Vb_w[0]) + (0.0 -
    -0.6 * open_loop_V1_B.u5)) + open_loop_V1_B.Mcg_b_idx_2 * open_loop_V1_B.Cn;
  open_loop_V1_B.c_theta_g[0] = open_loop_V1_B.c_theta_tmp_k;
  open_loop_V1_B.c_theta_g[1] = open_loop_V1_B.c_theta_tmp_c;
  open_loop_V1_B.c_theta_g[2] = open_loop_V1_B.c_theta_tmp_cx;
  open_loop_V1_B.c_theta_g[3] = open_loop_V1_B.c_theta_tmp_b;
  open_loop_V1_B.c_theta_g[4] = open_loop_V1_B.c_theta_tmp;
  open_loop_V1_B.c_theta_g[5] = open_loop_V1_B.c_theta_tmp_p;
  open_loop_V1_B.c_theta_g[6] = -open_loop_V1_B.s_theta;
  open_loop_V1_B.c_theta_g[7] = open_loop_V1_B.c_theta_tmp_cv;
  open_loop_V1_B.c_theta_g[8] = open_loop_V1_B.c_theta_tmp_f;
  open_loop_V1_B.CD_iw_IGE_m[0] = open_loop_V1_B.x[0];
  open_loop_V1_B.CD_iw_IGE_m[1] = open_loop_V1_B.x[1];
  open_loop_V1_B.CD_iw_IGE_m[2] = open_loop_V1_B.x[2];
  open_loop_V1_B.s_theta = open_loop_V1_B.Vb_w[0] + open_loop_V1_B.u5;
  open_loop_V1_B.F_b[0] = (open_loop_V1_B.Va + open_loop_V1_B.s_theta) +
    open_loop_V1_B.FA_b_idx_0;
  open_loop_V1_B.u5 = 0.0;
  open_loop_V1_B.F_b[1] = open_loop_V1_B.s_phi + open_loop_V1_B.FA_b_idx_1;
  open_loop_V1_B.F_b[2] = (open_loop_V1_B.c_phi + 0.17431148549531633) +
    open_loop_V1_B.FA_b_idx_2;
  open_loop_V1_B.Va_b_idx_2 = 0.0;
  open_loop_V1_B.Vb_w[0] = 39.71 * open_loop_V1_B.x[3] + 8.97 *
    open_loop_V1_B.x[5];
  open_loop_V1_B.Vb_w[1] = 85.51 * open_loop_V1_B.x[4];
  open_loop_V1_B.Vb_w[2] = 8.97 * open_loop_V1_B.x[3] + 114.39 *
    open_loop_V1_B.x[5];
  for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 3; open_loop_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&open_loop_V1_B.c_theta_g[3 *
      open_loop_V1_B.i]), _mm_set1_pd
      (open_loop_V1_B.CD_iw_IGE_m[open_loop_V1_B.i])), _mm_set_pd
                     (open_loop_V1_B.u4, open_loop_V1_B.u5));
    _mm_storeu_pd(&open_loop_V1_B.dv[0], tmp);
    open_loop_V1_B.u5 = open_loop_V1_B.dv[0];
    open_loop_V1_B.u4 = open_loop_V1_B.dv[1];
    open_loop_V1_B.Va_b_idx_2 += open_loop_V1_B.c_theta_g[3 * open_loop_V1_B.i +
      2] * open_loop_V1_B.CD_iw_IGE_m[open_loop_V1_B.i];
  }

  tmp = _mm_set_pd(open_loop_V1_B.x[3], open_loop_V1_B.x[5]);
  tmp_2 = _mm_sub_pd(_mm_mul_pd(_mm_set_pd(open_loop_V1_B.x[0],
    open_loop_V1_B.x[2]), _mm_loadu_pd(&open_loop_V1_B.x[4])), _mm_mul_pd
                     (_mm_loadu_pd(&open_loop_V1_B.x[1]), tmp));
  _mm_storeu_pd(&open_loop_V1_B.CD_iw_IGE_m[0], tmp_2);
  open_loop_V1_B.CD_iw_IGE_m[2] = open_loop_V1_B.x[1] * open_loop_V1_B.x[3] -
    open_loop_V1_B.x[0] * open_loop_V1_B.x[4];
  tmp = _mm_sub_pd(_mm_set_pd(open_loop_V1_B.Mcg_b_idx_1, open_loop_V1_B.c_theta),
                   _mm_sub_pd(_mm_mul_pd(_mm_set_pd(open_loop_V1_B.Vb_w[0],
    open_loop_V1_B.Vb_w[2]), _mm_loadu_pd(&open_loop_V1_B.x[4])), _mm_mul_pd
    (_mm_loadu_pd(&open_loop_V1_B.Vb_w[1]), tmp)));
  _mm_storeu_pd(&open_loop_V1_B.Mcg_b[0], tmp);
  open_loop_V1_B.Mcg_b[2] = open_loop_V1_B.Mcg_b_idx_2 - (open_loop_V1_B.Vb_w[1]
    * open_loop_V1_B.x[3] - open_loop_V1_B.Vb_w[0] * open_loop_V1_B.x[4]);
  open_loop_V1_B.c_theta_tmp_f = 0.0;
  open_loop_V1_B.c_theta_tmp_cv = 0.0;
  open_loop_V1_B.c_theta_tmp = 0.0;
  for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 3; open_loop_V1_B.i++) {
    _mm_storeu_pd(&open_loop_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&a[3
      * open_loop_V1_B.i]), _mm_set1_pd(open_loop_V1_B.Mcg_b[open_loop_V1_B.i])),
      _mm_set_pd(open_loop_V1_B.c_theta_tmp_cv, open_loop_V1_B.c_theta_tmp_f)));
    open_loop_V1_B.c_theta_tmp_f = open_loop_V1_B.dv[0];
    open_loop_V1_B.c_theta_tmp_cv = open_loop_V1_B.dv[1];
    open_loop_V1_B.c_theta_tmp += a_0[3 * open_loop_V1_B.i + 2] *
      open_loop_V1_B.Mcg_b[open_loop_V1_B.i];
    open_loop_V1_B.c_theta_g[3 * open_loop_V1_B.i] = c[open_loop_V1_B.i];
  }

  open_loop_V1_B.Vb_w[2] = open_loop_V1_B.c_theta_tmp;
  open_loop_V1_B.Vb_w[1] = open_loop_V1_B.c_theta_tmp_cv;
  open_loop_V1_B.Vb_w[0] = open_loop_V1_B.c_theta_tmp_f;
  open_loop_V1_B.c_theta_g[1] = 0.0;
  open_loop_V1_B.c_theta_g[4] = open_loop_V1_B.FA_b_tmp_m;
  open_loop_V1_B.c_theta_g[7] = -sin(open_loop_V1_B.s_psi);
  open_loop_V1_B.c_theta_g[2] = 0.0;
  open_loop_V1_B.c_theta_g[5] = open_loop_V1_B.FA_b_tmp;
  open_loop_V1_B.c_theta_g[8] = open_loop_V1_B.FA_b_tmp_m;
  open_loop_V1_B.c_theta_tmp_f = 0.0;
  open_loop_V1_B.c_theta_tmp_cv = 0.0;
  open_loop_V1_B.c_theta_tmp = 0.0;
  for (open_loop_V1_B.i = 0; open_loop_V1_B.i < 3; open_loop_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&open_loop_V1_B.c_theta_g[3 *
      open_loop_V1_B.i]), _mm_set1_pd(open_loop_V1_B.wbe_b[open_loop_V1_B.i])),
                     _mm_set_pd(open_loop_V1_B.c_theta_tmp_cv,
      open_loop_V1_B.c_theta_tmp_f));
    _mm_storeu_pd(&open_loop_V1_B.dv[0], tmp);
    open_loop_V1_B.c_theta_tmp_f = open_loop_V1_B.dv[0];
    open_loop_V1_B.c_theta_tmp_cv = open_loop_V1_B.dv[1];
    _mm_storeu_pd(&open_loop_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (0.0089285714285714281, open_loop_V1_B.c_theta_g[3 * open_loop_V1_B.i + 2]),
      _mm_set_pd(open_loop_V1_B.F_b[open_loop_V1_B.i],
                 open_loop_V1_B.wbe_b[open_loop_V1_B.i])), _mm_mul_pd(_mm_set_pd
      (open_loop_V1_B.CD_iw_IGE_m[open_loop_V1_B.i], open_loop_V1_B.c_theta_tmp),
      _mm_set_pd(-1.0, 1.0))));
    open_loop_V1_B.c_theta_tmp = open_loop_V1_B.dv[0];
    open_loop_V1_B.XDOT[open_loop_V1_B.i] = open_loop_V1_B.dv[1];
    open_loop_V1_B.XDOT[open_loop_V1_B.i + 3] =
      open_loop_V1_B.Vb_w[open_loop_V1_B.i];
  }

  open_loop_V1_B.XDOT[9] = open_loop_V1_B.u5;
  open_loop_V1_B.XDOT[10] = open_loop_V1_B.u4;
  open_loop_V1_B.XDOT[11] = open_loop_V1_B.Va_b_idx_2;
  open_loop_V1_B.XDOT[12] = open_loop_V1_B.CL_w_IGE / (open_loop_V1_B.CD_iw_IGE
    + 0.0306);
  open_loop_V1_B.XDOT[19] = 0.0;
  open_loop_V1_B.XDOT[20] = open_loop_V1_B.Cl;
  open_loop_V1_B.XDOT[21] = open_loop_V1_B.Cm;
  open_loop_V1_B.XDOT[22] = open_loop_V1_B.Cn;
  open_loop_V1_B.XDOT[23] = open_loop_V1_B.c_psi;
  open_loop_V1_B.XDOT[24] = open_loop_V1_B.s_psi;
  open_loop_V1_B.XDOT[25] = open_loop_V1_B.CL_w_OGE;
  open_loop_V1_B.XDOT[26] = open_loop_V1_B.CL_h_OGE;
  open_loop_V1_B.XDOT[27] = open_loop_V1_B.CL_w_IGE;
  open_loop_V1_B.XDOT[28] = open_loop_V1_B.u2;
  open_loop_V1_B.XDOT[29] = open_loop_V1_B.CD_iw_IGE;
  open_loop_V1_B.XDOT[30] = open_loop_V1_B.CD_ih_IGE;
  open_loop_V1_B.XDOT[6] = open_loop_V1_B.c_theta_tmp_f;
  open_loop_V1_B.XDOT[13] = open_loop_V1_B.F_b[0];
  open_loop_V1_B.XDOT[16] = open_loop_V1_B.c_theta;
  open_loop_V1_B.XDOT[31] = open_loop_V1_B.Va;
  open_loop_V1_B.XDOT[34] = open_loop_V1_B.s_theta;
  open_loop_V1_B.XDOT[37] = open_loop_V1_B.FA_b_idx_0;
  open_loop_V1_B.XDOT[7] = open_loop_V1_B.c_theta_tmp_cv;
  open_loop_V1_B.XDOT[14] = open_loop_V1_B.F_b[1];
  open_loop_V1_B.XDOT[17] = open_loop_V1_B.Mcg_b_idx_1;
  open_loop_V1_B.XDOT[32] = open_loop_V1_B.s_phi;
  open_loop_V1_B.XDOT[35] = 0.0;
  open_loop_V1_B.XDOT[38] = open_loop_V1_B.FA_b_idx_1;
  open_loop_V1_B.XDOT[8] = open_loop_V1_B.c_theta_tmp;
  open_loop_V1_B.XDOT[15] = open_loop_V1_B.F_b[2];
  open_loop_V1_B.XDOT[18] = open_loop_V1_B.Mcg_b_idx_2;
  open_loop_V1_B.XDOT[33] = open_loop_V1_B.c_phi;
  open_loop_V1_B.XDOT[36] = 0.17431148549531633;
  open_loop_V1_B.XDOT[39] = open_loop_V1_B.FA_b_idx_2;
  if (open_loop_V1_B.b) {
  }

  /* Gain: '<Root>/Gain' */
  open_loop_V1_B.Gain = -open_loop_V1_B.x[11];
  if (open_loop_V1_B.b) {
  }

  /* BusAssignment: '<Root>/Bus Assignment' */
  memset(&open_loop_V1_B.BusAssignment, 0, sizeof
         (SL_Bus_gazebo_msgs_SetEntityStateRequest));

  /* Gain: '<Root>/Gain2' incorporates:
   *  Gain: '<Root>/Gain1'
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   */
  _mm_storeu_pd(&open_loop_V1_B.dv[0], _mm_div_pd(_mm_set_pd(-open_loop_V1_B.x[7],
    -open_loop_V1_B.x[8]), _mm_set1_pd(2.0)));

  /* MATLABSystem: '<Root>/Coordinate Transformation Conversion' incorporates:
   *  Constant: '<Root>/Constant'
   *  Sum: '<Root>/Sum'
   */
  open_loop_V1_B.u2 = (open_loop_V1_B.x[6] + 1.5707963267948966) / 2.0;
  open_loop_V1_B.c_psi = sin(open_loop_V1_B.dv[0]);
  open_loop_V1_B.s_psi = sin(open_loop_V1_B.dv[1]);
  open_loop_V1_B.CL_w_OGE = sin(open_loop_V1_B.u2);
  open_loop_V1_B.CL_h_OGE = cos(open_loop_V1_B.dv[0]);
  open_loop_V1_B.CL_w_IGE = cos(open_loop_V1_B.dv[1]);
  open_loop_V1_B.u2 = cos(open_loop_V1_B.u2);

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  Gain: '<Root>/Gain3'
   */
  open_loop_V1_B.BusAssignment.state.pose.position.x = open_loop_V1_B.x[9];
  open_loop_V1_B.BusAssignment.state.pose.position.y = -open_loop_V1_B.x[10];
  open_loop_V1_B.BusAssignment.state.pose.position.z = open_loop_V1_B.Gain;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  open_loop_V1_B.CD_iw_IGE = open_loop_V1_B.CL_h_OGE * open_loop_V1_B.CL_w_IGE;

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   * */
  open_loop_V1_B.BusAssignment.state.pose.orientation.w = open_loop_V1_B.c_psi *
    open_loop_V1_B.s_psi * open_loop_V1_B.CL_w_OGE + open_loop_V1_B.CD_iw_IGE *
    open_loop_V1_B.u2;
  open_loop_V1_B.BusAssignment.state.pose.orientation.z =
    open_loop_V1_B.CD_iw_IGE * open_loop_V1_B.CL_w_OGE - open_loop_V1_B.u2 *
    open_loop_V1_B.c_psi * open_loop_V1_B.s_psi;
  open_loop_V1_B.BusAssignment.state.pose.orientation.y =
    open_loop_V1_B.CL_h_OGE * open_loop_V1_B.u2 * open_loop_V1_B.s_psi +
    open_loop_V1_B.CL_w_IGE * open_loop_V1_B.c_psi * open_loop_V1_B.CL_w_OGE;
  open_loop_V1_B.BusAssignment.state.pose.orientation.x =
    open_loop_V1_B.CL_w_IGE * open_loop_V1_B.u2 * open_loop_V1_B.c_psi -
    open_loop_V1_B.CL_h_OGE * open_loop_V1_B.s_psi * open_loop_V1_B.CL_w_OGE;
  memcpy(&open_loop_V1_B.BusAssignment.state.name[0],
         &open_loop_V1_B.stringOut_l[0], sizeof(uint8_T) << 7U);
  memcpy(&open_loop_V1_B.BusAssignment.state.reference_frame[0],
         &open_loop_V1_B.stringOut[0], sizeof(uint8_T) << 7U);
  open_loop_V1_B.BusAssignment.state.name_SL_Info.CurrentLength =
    open_loop_V1_B.lengthOut_e;
  open_loop_V1_B.BusAssignment.state.reference_frame_SL_Info.CurrentLength =
    open_loop_V1_B.lengthOut;

  /* Outputs for Atomic SubSystem: '<Root>/Call Service' */
  /* MATLABSystem: '<S2>/ServiceCaller' */
  open_loop_V1_B.serverAvailableOnTime = ServCall_open_loop_V1_326.waitForServer
    (5.0);
  if (open_loop_V1_B.serverAvailableOnTime) {
    ServCall_open_loop_V1_326.call(&open_loop_V1_B.BusAssignment,
      &open_loop_V1_B.r);
  }

  /* End of MATLABSystem: '<S2>/ServiceCaller' */
  /* End of Outputs for SubSystem: '<Root>/Call Service' */

  /* Product: '<S7>/Product2' incorporates:
   *  Math: '<S7>/Square'
   *  Math: '<S7>/Square1'
   *  Math: '<S7>/Square2'
   *  Sqrt: '<S7>/Sqrt'
   *  Sum: '<S7>/Sum2'
   */
  open_loop_V1_B.Power = sqrt((open_loop_V1_B.x[0] * open_loop_V1_B.x[0] +
    open_loop_V1_B.x[1] * open_loop_V1_B.x[1]) + open_loop_V1_B.x[2] *
    open_loop_V1_B.x[2]) * open_loop_V1_B.XDOT[34];

  /* Gain: '<S7>/Gain3' */
  open_loop_V1_B.Gain3 = 0.001 * open_loop_V1_B.Power;
  if (open_loop_V1_B.b) {
  }

  /* Product: '<S7>/Divide' incorporates:
   *  Constant: '<S7>/thrust efficiency Cp?'
   */
  open_loop_V1_B.powerdemand = open_loop_V1_B.Gain3 / 0.248;
  if (open_loop_V1_B.b) {
  }

  /* Product: '<S7>/Divide1' */
  open_loop_V1_B.loadtorque = open_loop_V1_B.powerdemand /
    open_loop_V1_ConstB.motorspeed;
  if (open_loop_V1_B.b) {
  }

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Constant: '<Root>/Constant2'
   *  Gain: '<S39>/Integral Gain'
   *  Gain: '<S91>/Integral Gain'
   *  Sum: '<Root>/Sum1'
   */
  _mm_storeu_pd(&open_loop_V1_B.dv[0], _mm_mul_pd(_mm_sub_pd(_mm_set_pd(1.0,
    -0.034906585039886591), _mm_set_pd(open_loop_V1_B.x[11], open_loop_V1_B.x[7])),
    _mm_set_pd(1.0E-6, 0.01)));

  /* Gain: '<S39>/Integral Gain' */
  open_loop_V1_B.IntegralGain = open_loop_V1_B.dv[0];

  /* Gain: '<S91>/Integral Gain' */
  open_loop_V1_B.IntegralGain_m = open_loop_V1_B.dv[1];

  /* Gain: '<S7>/Gain1' incorporates:
   *  Integrator: '<S7>/Integrator1'
   */
  open_loop_V1_B.EnergykWh = 2.7777777777777776E-7 *
    open_loop_V1_X.Integrator1_CSTATE;
  if (open_loop_V1_B.b) {
  }

  if (rtmIsMajorTimeStep((&open_loop_V1_M))) {
    rt_ertODEUpdateContinuousStates(&(&open_loop_V1_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&open_loop_V1_M)->Timing.clockTick0)) {
      ++(&open_loop_V1_M)->Timing.clockTickH0;
    }

    (&open_loop_V1_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&open_loop_V1_M)
      ->solverInfo);

    {
      /* Update absolute timer for sample time: [0.01s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.01, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      (&open_loop_V1_M)->Timing.clockTick1++;
      if (!(&open_loop_V1_M)->Timing.clockTick1) {
        (&open_loop_V1_M)->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void open_loop_V1::open_loop_V1_derivatives()
{
  XDot_open_loop_V1_T *_rtXdot;
  _rtXdot = ((XDot_open_loop_V1_T *) (&open_loop_V1_M)->derivs);

  /* Derivatives for Integrator: '<S7>/Integrator' */
  memcpy(&_rtXdot->Integrator_CSTATE[0], &open_loop_V1_B.XDOT[0], 12U * sizeof
         (real_T));

  /* Derivatives for Integrator: '<S94>/Integrator' */
  _rtXdot->Integrator_CSTATE_b = open_loop_V1_B.IntegralGain_m;

  /* Derivatives for Integrator: '<S89>/Filter' */
  _rtXdot->Filter_CSTATE = open_loop_V1_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S42>/Integrator' */
  _rtXdot->Integrator_CSTATE_p = open_loop_V1_B.IntegralGain;

  /* Derivatives for Integrator: '<S37>/Filter' */
  _rtXdot->Filter_CSTATE_m = open_loop_V1_B.FilterCoefficient_m;

  /* Derivatives for Integrator: '<S7>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = open_loop_V1_B.Power;
}

/* Model initialize function */
void open_loop_V1::initialize()
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&open_loop_V1_M)->solverInfo, &(&open_loop_V1_M)
                          ->Timing.simTimeStep);
    rtsiSetTPtr(&(&open_loop_V1_M)->solverInfo, &rtmGetTPtr((&open_loop_V1_M)));
    rtsiSetStepSizePtr(&(&open_loop_V1_M)->solverInfo, &(&open_loop_V1_M)
                       ->Timing.stepSize0);
    rtsiSetdXPtr(&(&open_loop_V1_M)->solverInfo, &(&open_loop_V1_M)->derivs);
    rtsiSetContStatesPtr(&(&open_loop_V1_M)->solverInfo, (real_T **)
                         &(&open_loop_V1_M)->contStates);
    rtsiSetNumContStatesPtr(&(&open_loop_V1_M)->solverInfo, &(&open_loop_V1_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&open_loop_V1_M)->solverInfo,
      &(&open_loop_V1_M)->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&open_loop_V1_M)->solverInfo,
      &(&open_loop_V1_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&open_loop_V1_M)->solverInfo,
      &(&open_loop_V1_M)->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&open_loop_V1_M)->solverInfo, (boolean_T**) &(
      &open_loop_V1_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&open_loop_V1_M)->solverInfo, (&rtmGetErrorStatus
      ((&open_loop_V1_M))));
    rtsiSetRTModelPtr(&(&open_loop_V1_M)->solverInfo, (&open_loop_V1_M));
  }

  rtsiSetSimTimeStep(&(&open_loop_V1_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&open_loop_V1_M)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&open_loop_V1_M)->solverInfo, false);
  (&open_loop_V1_M)->intgData.y = (&open_loop_V1_M)->odeY;
  (&open_loop_V1_M)->intgData.f[0] = (&open_loop_V1_M)->odeF[0];
  (&open_loop_V1_M)->intgData.f[1] = (&open_loop_V1_M)->odeF[1];
  (&open_loop_V1_M)->intgData.f[2] = (&open_loop_V1_M)->odeF[2];
  (&open_loop_V1_M)->contStates = ((X_open_loop_V1_T *) &open_loop_V1_X);
  (&open_loop_V1_M)->contStateDisabled = ((XDis_open_loop_V1_T *)
    &open_loop_V1_XDis);
  (&open_loop_V1_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&open_loop_V1_M)->solverInfo, static_cast<void *>
                    (&(&open_loop_V1_M)->intgData));
  rtsiSetSolverName(&(&open_loop_V1_M)->solverInfo,"ode3");
  rtmSetTPtr((&open_loop_V1_M), &(&open_loop_V1_M)->Timing.tArray[0]);
  (&open_loop_V1_M)->Timing.stepSize0 = 0.01;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  open_loop_V1_DW.objisempty = true;
  open_loop_V1_DW.obj_c.isInitialized = 1;

  /* Start for Atomic SubSystem: '<Root>/Call Service' */
  /* Start for MATLABSystem: '<S2>/ServiceCaller' */
  open_loop_V1_DW.obj.QOSAvoidROSNamespaceConventions = false;
  open_loop_V1_DW.obj.matlabCodegenIsDeleted = false;
  open_loop_V1_DW.objisempty_f = true;
  open_loop_V1_DW.obj.isSetupComplete = false;
  open_loop_V1_DW.obj.isInitialized = 1;
  open_lo_ServiceCaller_setupImpl(&open_loop_V1_DW.obj);
  open_loop_V1_DW.obj.isSetupComplete = true;

  /* End of Start for SubSystem: '<Root>/Call Service' */

  /* InitializeConditions for Integrator: '<S7>/Integrator' */
  memcpy(&open_loop_V1_X.Integrator_CSTATE[0],
         &open_loop_V1_ConstP.Integrator_IC[0], 12U * sizeof(real_T));

  /* InitializeConditions for Integrator: '<S94>/Integrator' */
  open_loop_V1_X.Integrator_CSTATE_b = 0.0;

  /* InitializeConditions for Integrator: '<S89>/Filter' */
  open_loop_V1_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S42>/Integrator' */
  open_loop_V1_X.Integrator_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S37>/Filter' */
  open_loop_V1_X.Filter_CSTATE_m = 0.0;

  /* InitializeConditions for Integrator: '<S7>/Integrator1' */
  open_loop_V1_X.Integrator1_CSTATE = 0.0;
}

/* Model terminate function */
void open_loop_V1::terminate()
{
  /* Terminate for Atomic SubSystem: '<Root>/Call Service' */
  /* Terminate for MATLABSystem: '<S2>/ServiceCaller' */
  if (!open_loop_V1_DW.obj.matlabCodegenIsDeleted) {
    open_loop_V1_DW.obj.matlabCodegenIsDeleted = true;
    if ((open_loop_V1_DW.obj.isInitialized == 1) &&
        open_loop_V1_DW.obj.isSetupComplete) {
      ServCall_open_loop_V1_326.resetSvcClientPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S2>/ServiceCaller' */
  /* End of Terminate for SubSystem: '<Root>/Call Service' */
}

/* Constructor */
open_loop_V1::open_loop_V1() :
  open_loop_V1_B(),
  open_loop_V1_DW(),
  open_loop_V1_X(),
  open_loop_V1_XDis(),
  open_loop_V1_M()
{
  /* Currently there is no constructor body generated.*/
}

/* Destructor */
open_loop_V1::~open_loop_V1()
{
  /* Currently there is no destructor body generated.*/
}

/* Real-Time Model get method */
RT_MODEL_open_loop_V1_T * open_loop_V1::getRTM()
{
  return (&open_loop_V1_M);
}
