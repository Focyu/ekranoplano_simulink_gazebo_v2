/*
 * pid_control_V1.cpp
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.27
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Fri Feb 27 22:52:44 2026
 *
 * Target selection: ert.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "pid_control_V1.h"
#include "pid_control_V1_types.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include <emmintrin.h>

extern "C"
{

#include "rt_nonfinite.h"

}

#include "rmw/qos_profiles.h"
#include <stddef.h>
#include "rt_defines.h"
#include "pid_control_V1_private.h"

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
void pid_control_V1::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
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
  pid_control_V1_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  this->step();
  pid_control_V1_derivatives();

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
  pid_control_V1_derivatives();

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

void pid_control_V1::pid_con_ServiceCaller_setupImpl(const
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
    pid_control_V1_B.b_zeroDelimTopic[i] = b_zeroDelimTopic[i];
  }

  ServCall_pid_control_V1_326.createServiceCaller
    (&pid_control_V1_B.b_zeroDelimTopic[0], qos_profile);
}

real_T pid_control_V1::pid_control_V1_rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      pid_control_V1_B.i1 = 1;
    } else {
      pid_control_V1_B.i1 = -1;
    }

    if (u1 > 0.0) {
      pid_control_V1_B.i2 = 1;
    } else {
      pid_control_V1_B.i2 = -1;
    }

    y = atan2(static_cast<real_T>(pid_control_V1_B.i1), static_cast<real_T>
              (pid_control_V1_B.i2));
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

real_T pid_control_V1::pid_control_V1_rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    pid_control_V1_B.d = fabs(u0);
    pid_control_V1_B.d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (pid_control_V1_B.d == 1.0) {
        y = 1.0;
      } else if (pid_control_V1_B.d > 1.0) {
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
    } else if (pid_control_V1_B.d1 == 0.0) {
      y = 1.0;
    } else if (pid_control_V1_B.d1 == 1.0) {
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
void pid_control_V1::step()
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

  if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
    /* set solver stop time */
    if (!((&pid_control_V1_M)->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&(&pid_control_V1_M)->solverInfo,
                            (((&pid_control_V1_M)->Timing.clockTickH0 + 1) *
        (&pid_control_V1_M)->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&(&pid_control_V1_M)->solverInfo,
                            (((&pid_control_V1_M)->Timing.clockTick0 + 1) *
        (&pid_control_V1_M)->Timing.stepSize0 + (&pid_control_V1_M)
        ->Timing.clockTickH0 * (&pid_control_V1_M)->Timing.stepSize0 *
        4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&pid_control_V1_M))) {
    (&pid_control_V1_M)->Timing.t[0] = rtsiGetT(&(&pid_control_V1_M)->solverInfo);
  }

  pid_control_V1_B.b = rtmIsMajorTimeStep((&pid_control_V1_M));
  if (pid_control_V1_B.b) {
    /* MATLAB Function: '<Root>/MATLAB Function' */
    memset(&pid_control_V1_B.stringOut_l[0], 0, sizeof(uint8_T) << 7U);
    for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 11; pid_control_V1_B.i++)
    {
      pid_control_V1_B.stringOut_l[pid_control_V1_B.i] = b[pid_control_V1_B.i];
    }

    pid_control_V1_B.lengthOut_e = 11U;

    /* End of MATLAB Function: '<Root>/MATLAB Function' */

    /* MATLAB Function: '<Root>/MATLAB Function1' */
    memset(&pid_control_V1_B.stringOut[0], 0, sizeof(uint8_T) << 7U);
    for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 5; pid_control_V1_B.i++) {
      pid_control_V1_B.stringOut[pid_control_V1_B.i] = b_0[pid_control_V1_B.i];
    }

    pid_control_V1_B.lengthOut = 5U;

    /* End of MATLAB Function: '<Root>/MATLAB Function1' */
  }

  /* Integrator: '<S7>/Integrator' */
  memcpy(&pid_control_V1_B.x[0], &pid_control_V1_X.Integrator_CSTATE[0], 12U *
         sizeof(real_T));

  /* Gain: '<S47>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S37>/Derivative Gain'
   *  Integrator: '<S39>/Filter'
   *  Sum: '<Root>/Sum1'
   *  Sum: '<S39>/SumD'
   */
  pid_control_V1_B.FilterCoefficient = ((0.0087266462599716477 -
    pid_control_V1_B.x[7]) * 4.5 - pid_control_V1_X.Filter_CSTATE) * 100.0;

  /* Sum: '<S53>/Sum' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S49>/Proportional Gain'
   *  Integrator: '<S44>/Integrator'
   *  Sum: '<Root>/Sum1'
   */
  pid_control_V1_B.SignPreSat = ((0.0087266462599716477 - pid_control_V1_B.x[7])
    * 8.0 + pid_control_V1_X.Integrator_CSTATE_p) +
    pid_control_V1_B.FilterCoefficient;

  /* Saturate: '<S51>/Saturation' */
  if (pid_control_V1_B.SignPreSat > 0.3490658503988659) {
    /* Saturate: '<S51>/Saturation' */
    pid_control_V1_B.Saturation = 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat < -0.3490658503988659) {
    /* Saturate: '<S51>/Saturation' */
    pid_control_V1_B.Saturation = -0.3490658503988659;
  } else {
    /* Saturate: '<S51>/Saturation' */
    pid_control_V1_B.Saturation = pid_control_V1_B.SignPreSat;
  }

  /* End of Saturate: '<S51>/Saturation' */

  /* Gain: '<S101>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Gain: '<Root>/Gain4'
   *  Gain: '<S91>/Derivative Gain'
   *  Integrator: '<S93>/Filter'
   *  Sum: '<Root>/Sum2'
   *  Sum: '<S93>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_g = ((1.0 - (-pid_control_V1_B.x[11])) *
    0.15 - pid_control_V1_X.Filter_CSTATE_f) * 100.0;

  /* Sum: '<S107>/Sum' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Gain: '<Root>/Gain4'
   *  Gain: '<S103>/Proportional Gain'
   *  Integrator: '<S98>/Integrator'
   *  Sum: '<Root>/Sum2'
   */
  pid_control_V1_B.SignPreSat_o = ((1.0 - (-pid_control_V1_B.x[11])) * 0.1 +
    pid_control_V1_X.Integrator_CSTATE_b) + pid_control_V1_B.FilterCoefficient_g;

  /* Saturate: '<S105>/Saturation' */
  if (pid_control_V1_B.SignPreSat_o > 1.0) {
    /* Saturate: '<S105>/Saturation' */
    pid_control_V1_B.Saturation_p = 1.0;
  } else if (pid_control_V1_B.SignPreSat_o < 0.0) {
    /* Saturate: '<S105>/Saturation' */
    pid_control_V1_B.Saturation_p = 0.0;
  } else {
    /* Saturate: '<S105>/Saturation' */
    pid_control_V1_B.Saturation_p = pid_control_V1_B.SignPreSat_o;
  }

  /* End of Saturate: '<S105>/Saturation' */

  /* Gain: '<Root>/Gain5' */
  pid_control_V1_B.Gain5 = 0.000125 * pid_control_V1_B.Saturation_p;

  /* SignalConversion generated from: '<S116>/ SFunction ' incorporates:
   *  Constant: '<Root>/nominal_control'
   *  MATLAB Function: '<S7>/MATLAB Function'
   */
  pid_control_V1_B.control_vector[0] = 0.0;
  pid_control_V1_B.control_vector[1] = pid_control_V1_B.Saturation;
  pid_control_V1_B.control_vector[2] = 0.0;
  pid_control_V1_B.control_vector[3] = pid_control_V1_B.Gain5;
  pid_control_V1_B.control_vector[4] = pid_control_V1_B.Gain5;

  /* MATLAB Function: '<S7>/MATLAB Function' */
  pid_control_V1_B.u2 = pid_control_V1_B.control_vector[1];
  pid_control_V1_B.u4 = pid_control_V1_B.control_vector[3];
  pid_control_V1_B.u5 = pid_control_V1_B.control_vector[4];
  if (pid_control_V1_B.control_vector[1] > 0.3490658503988659) {
    pid_control_V1_B.u2 = 0.3490658503988659;
  } else if (pid_control_V1_B.control_vector[1] < -0.3490658503988659) {
    pid_control_V1_B.u2 = -0.3490658503988659;
  }

  if (pid_control_V1_B.control_vector[3] > 0.000125) {
    pid_control_V1_B.u4 = 0.000125;
  } else if (pid_control_V1_B.control_vector[3] < 0.0) {
    pid_control_V1_B.u4 = 0.0;
  }

  if (pid_control_V1_B.control_vector[4] > 0.000125) {
    pid_control_V1_B.u5 = 0.000125;
  } else if (pid_control_V1_B.control_vector[4] < 0.0) {
    pid_control_V1_B.u5 = 0.0;
  }

  pid_control_V1_B.c_phi = cos(pid_control_V1_B.x[6]);
  pid_control_V1_B.s_phi = sin(pid_control_V1_B.x[6]);
  pid_control_V1_B.c_theta_m = cos(pid_control_V1_B.x[7]);
  pid_control_V1_B.s_theta = sin(pid_control_V1_B.x[7]);
  pid_control_V1_B.c_psi = cos(pid_control_V1_B.x[8]);
  pid_control_V1_B.s_psi = sin(pid_control_V1_B.x[8]);
  pid_control_V1_B.c_theta_tmp_g1 = pid_control_V1_B.c_theta_m *
    pid_control_V1_B.c_psi;
  pid_control_V1_B.c_theta[0] = pid_control_V1_B.c_theta_tmp_g1;
  pid_control_V1_B.c_theta_tmp_g = pid_control_V1_B.c_theta_m *
    pid_control_V1_B.s_psi;
  pid_control_V1_B.c_theta[3] = pid_control_V1_B.c_theta_tmp_g;
  pid_control_V1_B.c_theta[6] = -pid_control_V1_B.s_theta;
  pid_control_V1_B.c_theta_tmp = pid_control_V1_B.s_phi *
    pid_control_V1_B.s_theta;
  pid_control_V1_B.c_theta_tmp_f = pid_control_V1_B.c_theta_tmp *
    pid_control_V1_B.c_psi - pid_control_V1_B.c_phi * pid_control_V1_B.s_psi;
  pid_control_V1_B.c_theta[1] = pid_control_V1_B.c_theta_tmp_f;
  pid_control_V1_B.c_theta_tmp = pid_control_V1_B.c_theta_tmp *
    pid_control_V1_B.s_psi + pid_control_V1_B.c_phi * pid_control_V1_B.c_psi;
  pid_control_V1_B.c_theta[4] = pid_control_V1_B.c_theta_tmp;
  pid_control_V1_B.c_theta_tmp_cv = pid_control_V1_B.s_phi *
    pid_control_V1_B.c_theta_m;
  pid_control_V1_B.c_theta[7] = pid_control_V1_B.c_theta_tmp_cv;
  pid_control_V1_B.c_theta_tmp_c = pid_control_V1_B.c_phi *
    pid_control_V1_B.s_theta;
  pid_control_V1_B.c_theta_tmp_p = pid_control_V1_B.c_theta_tmp_c *
    pid_control_V1_B.c_psi + pid_control_V1_B.s_phi * pid_control_V1_B.s_psi;
  pid_control_V1_B.c_theta[2] = pid_control_V1_B.c_theta_tmp_p;
  pid_control_V1_B.c_theta_tmp_c = pid_control_V1_B.c_theta_tmp_c *
    pid_control_V1_B.s_psi - pid_control_V1_B.s_phi * pid_control_V1_B.c_psi;
  pid_control_V1_B.c_theta[5] = pid_control_V1_B.c_theta_tmp_c;
  pid_control_V1_B.c_theta_tmp_b = pid_control_V1_B.c_phi *
    pid_control_V1_B.c_theta_m;
  pid_control_V1_B.c_theta[8] = pid_control_V1_B.c_theta_tmp_b;
  for (pid_control_V1_B.i = 0; pid_control_V1_B.i <= 0; pid_control_V1_B.i += 2)
  {
    tmp = _mm_loadu_pd(&pid_control_V1_B.c_theta[pid_control_V1_B.i + 3]);
    tmp_2 = _mm_set1_pd(0.0);
    tmp_0 = _mm_loadu_pd(&pid_control_V1_B.c_theta[pid_control_V1_B.i]);
    tmp_1 = _mm_loadu_pd(&pid_control_V1_B.c_theta[pid_control_V1_B.i + 6]);
    _mm_storeu_pd(&pid_control_V1_B.Vb_w[pid_control_V1_B.i], _mm_add_pd
                  (_mm_add_pd(_mm_mul_pd(tmp, tmp_2), _mm_mul_pd(tmp_0, tmp_2)),
                   _mm_mul_pd(tmp_1, tmp_2)));
  }

  for (pid_control_V1_B.i = 2; pid_control_V1_B.i < 3; pid_control_V1_B.i++) {
    pid_control_V1_B.Vb_w[pid_control_V1_B.i] =
      (pid_control_V1_B.c_theta[pid_control_V1_B.i + 3] * 0.0 +
       pid_control_V1_B.c_theta[pid_control_V1_B.i] * 0.0) +
      pid_control_V1_B.c_theta[pid_control_V1_B.i + 6] * 0.0;
  }

  tmp = _mm_sub_pd(_mm_loadu_pd(&pid_control_V1_B.x[0]), _mm_loadu_pd
                   (&pid_control_V1_B.Vb_w[0]));
  _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp);

  /* MATLAB Function: '<S7>/MATLAB Function' */
  pid_control_V1_B.Va_b_idx_2 = pid_control_V1_B.x[2] - pid_control_V1_B.Vb_w[2];
  pid_control_V1_B.Va = sqrt((pid_control_V1_B.dv[0] * pid_control_V1_B.dv[0] +
    pid_control_V1_B.dv[1] * pid_control_V1_B.dv[1]) +
    pid_control_V1_B.Va_b_idx_2 * pid_control_V1_B.Va_b_idx_2);
  pid_control_V1_B.c_psi = pid_control_V1_rt_atan2d_snf
    (pid_control_V1_B.Va_b_idx_2, pid_control_V1_B.dv[0]);
  pid_control_V1_B.s_psi = asin(pid_control_V1_B.dv[1] / pid_control_V1_B.Va);
  pid_control_V1_B.Va_b_idx_2 = pid_control_V1_B.Va * pid_control_V1_B.Va *
    0.6125;
  pid_control_V1_B.wbe_b[0] = pid_control_V1_B.x[3];
  pid_control_V1_B.wbe_b[1] = pid_control_V1_B.x[4];
  pid_control_V1_B.wbe_b[2] = pid_control_V1_B.x[5];
  pid_control_V1_B.CL_w_OGE = ((pid_control_V1_B.c_psi - -0.065449846949787352)
    + 0.082903139469730644) * 4.9604094530365153;
  pid_control_V1_B.CL_h_OGE = ((pid_control_V1_B.c_psi - -0.074176493209759012)
    + 0.043633231299858237) * 4.8387748917360032;
  pid_control_V1_B.CD_iw_IGE_c = fabs((-pid_control_V1_B.x[11] - 0.363) / 5.02);
  pid_control_V1_B.CL_w_IGE = (288.0 * pid_control_V1_rt_powd_snf
    (pid_control_V1_B.CD_iw_IGE_c, 0.787) * exp(-9.14 *
    pid_control_V1_rt_powd_snf(pid_control_V1_B.CD_iw_IGE_c, 0.327)) *
    0.97986308862072491 / 5.9129476540958859 + 1.0) * pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.CD_ih_IGE = fabs((-pid_control_V1_B.x[11] + 0.72) / 2.74);
  pid_control_V1_B.u2 = (288.0 * pid_control_V1_rt_powd_snf
    (pid_control_V1_B.CD_ih_IGE, 0.787) * exp(-9.14 * pid_control_V1_rt_powd_snf
    (pid_control_V1_B.CD_ih_IGE, 0.327)) * 0.95628590200128227 /
    5.35300902982722 + 1.0) * ((pid_control_V1_B.u2 * pid_control_V1_B.u2 *
    -0.00141 + 0.0307 * pid_control_V1_B.u2) + pid_control_V1_B.CL_h_OGE);
  pid_control_V1_B.CD_iw_IGE_c = pid_control_V1_B.CL_w_OGE *
    pid_control_V1_B.CL_w_OGE / 21.205750411731103 * (1.0 - exp(-10.1 *
    pid_control_V1_rt_powd_snf(pid_control_V1_B.CD_iw_IGE_c, 0.686)));
  pid_control_V1_B.CD_ih_IGE = pid_control_V1_B.CL_h_OGE *
    pid_control_V1_B.CL_h_OGE / 18.943803701146454 * (1.0 - exp(-10.1 *
    pid_control_V1_rt_powd_snf(pid_control_V1_B.CD_ih_IGE, 0.686)));
  pid_control_V1_B.FA_b_tmp = sin(pid_control_V1_B.c_psi);
  pid_control_V1_B.FA_b_tmp_k = cos(pid_control_V1_B.c_psi);
  pid_control_V1_B.c_theta[0] = pid_control_V1_B.FA_b_tmp_k;
  pid_control_V1_B.c_theta[3] = 0.0;
  pid_control_V1_B.c_theta[6] = -pid_control_V1_B.FA_b_tmp;
  pid_control_V1_B.c_theta[2] = pid_control_V1_B.FA_b_tmp;
  pid_control_V1_B.c_theta[5] = 0.0;
  pid_control_V1_B.c_theta[8] = pid_control_V1_B.FA_b_tmp_k;
  pid_control_V1_B.CD_iw_IGE[0] = -(((pid_control_V1_B.CD_iw_IGE_c + 0.0306) +
    (pid_control_V1_B.CD_ih_IGE + 0.0306)) * pid_control_V1_B.Va_b_idx_2 * 3.334);
  pid_control_V1_B.CD_iw_IGE[1] = 0.0 * pid_control_V1_B.Va_b_idx_2 * 3.334;
  pid_control_V1_B.CD_iw_IGE[2] = -((pid_control_V1_B.CL_w_IGE * 3.334 +
    pid_control_V1_B.u2 * 1.128) * pid_control_V1_B.Va_b_idx_2);
  pid_control_V1_B.c_theta[1] = 0.0;
  pid_control_V1_B.FA_b_idx_0 = 0.0;
  pid_control_V1_B.c_theta[4] = 1.0;
  pid_control_V1_B.FA_b_idx_1 = 0.0;
  pid_control_V1_B.c_theta[7] = 0.0;
  pid_control_V1_B.FA_b_idx_2 = 0.0;
  for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 3; pid_control_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.c_theta[3 *
      pid_control_V1_B.i]), _mm_set1_pd
      (pid_control_V1_B.CD_iw_IGE[pid_control_V1_B.i])), _mm_set_pd
                     (pid_control_V1_B.FA_b_idx_1, pid_control_V1_B.FA_b_idx_0));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp);
    pid_control_V1_B.FA_b_idx_0 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FA_b_idx_1 = pid_control_V1_B.dv[1];
    pid_control_V1_B.FA_b_idx_2 += pid_control_V1_B.c_theta[3 *
      pid_control_V1_B.i + 2] * pid_control_V1_B.CD_iw_IGE[pid_control_V1_B.i];
  }

  pid_control_V1_B.Cl = -0.0005 * pid_control_V1_B.s_psi * 180.0 /
    3.1415926535897931;
  pid_control_V1_B.Cm = -0.02 * pid_control_V1_B.c_psi * 180.0 /
    3.1415926535897931;
  pid_control_V1_B.Cn = -0.002 * pid_control_V1_B.s_psi * 180.0 /
    3.1415926535897931;
  pid_control_V1_B.u4 = pid_control_V1_B.u4 / 0.000125 * (37.42 -
    pid_control_V1_B.Va) + pid_control_V1_B.Va;
  pid_control_V1_B.u5 = pid_control_V1_B.u5 / 0.000125 * (37.42 -
    pid_control_V1_B.Va) + pid_control_V1_B.Va;
  pid_control_V1_B.Vb_w[0] = 0.11056515539994617 * pid_control_V1_B.u4 *
    (pid_control_V1_B.u4 - pid_control_V1_B.Va) * 0.99619469809174555;
  pid_control_V1_B.u5 = 0.11056515539994617 * pid_control_V1_B.u5 *
    (pid_control_V1_B.u5 - pid_control_V1_B.Va) * 0.99619469809174555;
  pid_control_V1_B.u4 = 0.0;
  pid_control_V1_B.Va = -9.81 * pid_control_V1_B.s_theta * 112.0;
  pid_control_V1_B.c_theta_m *= 9.81;
  pid_control_V1_B.s_phi = pid_control_V1_B.c_theta_m * pid_control_V1_B.s_phi *
    112.0;
  pid_control_V1_B.c_phi = pid_control_V1_B.c_theta_m * pid_control_V1_B.c_phi *
    112.0;
  pid_control_V1_B.Mcg_b_idx_2 = 5.02 * pid_control_V1_B.Va_b_idx_2 * 3.334;
  pid_control_V1_B.c_theta_m = pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cl;
  pid_control_V1_B.Mcg_b_idx_1 = 0.646 * pid_control_V1_B.Va_b_idx_2 * 3.334 *
    pid_control_V1_B.Cm + ((-0.285 * pid_control_V1_B.Vb_w[0] -
    0.04523383048603459) + (-0.285 * pid_control_V1_B.u5 - 0.04523383048603459));
  pid_control_V1_B.Mcg_b_idx_2 = ((0.0 - 0.6 * pid_control_V1_B.Vb_w[0]) + (0.0
    - -0.6 * pid_control_V1_B.u5)) + pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cn;
  pid_control_V1_B.c_theta[0] = pid_control_V1_B.c_theta_tmp_b;
  pid_control_V1_B.c_theta[1] = pid_control_V1_B.c_theta_tmp_c;
  pid_control_V1_B.c_theta[2] = pid_control_V1_B.c_theta_tmp_p;
  pid_control_V1_B.c_theta[3] = pid_control_V1_B.c_theta_tmp_cv;
  pid_control_V1_B.c_theta[4] = pid_control_V1_B.c_theta_tmp;
  pid_control_V1_B.c_theta[5] = pid_control_V1_B.c_theta_tmp_f;
  pid_control_V1_B.c_theta[6] = -pid_control_V1_B.s_theta;
  pid_control_V1_B.c_theta[7] = pid_control_V1_B.c_theta_tmp_g;
  pid_control_V1_B.c_theta[8] = pid_control_V1_B.c_theta_tmp_g1;
  pid_control_V1_B.CD_iw_IGE[0] = pid_control_V1_B.x[0];
  pid_control_V1_B.CD_iw_IGE[1] = pid_control_V1_B.x[1];
  pid_control_V1_B.CD_iw_IGE[2] = pid_control_V1_B.x[2];
  pid_control_V1_B.s_theta = pid_control_V1_B.Vb_w[0] + pid_control_V1_B.u5;
  pid_control_V1_B.F_b[0] = (pid_control_V1_B.Va + pid_control_V1_B.s_theta) +
    pid_control_V1_B.FA_b_idx_0;
  pid_control_V1_B.u5 = 0.0;
  pid_control_V1_B.F_b[1] = pid_control_V1_B.s_phi + pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.F_b[2] = (pid_control_V1_B.c_phi + 0.17431148549531633) +
    pid_control_V1_B.FA_b_idx_2;
  pid_control_V1_B.Va_b_idx_2 = 0.0;
  pid_control_V1_B.Vb_w[0] = 39.71 * pid_control_V1_B.x[3] + 8.97 *
    pid_control_V1_B.x[5];
  pid_control_V1_B.Vb_w[1] = 85.51 * pid_control_V1_B.x[4];
  pid_control_V1_B.Vb_w[2] = 8.97 * pid_control_V1_B.x[3] + 114.39 *
    pid_control_V1_B.x[5];
  for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 3; pid_control_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.c_theta[3 *
      pid_control_V1_B.i]), _mm_set1_pd
      (pid_control_V1_B.CD_iw_IGE[pid_control_V1_B.i])), _mm_set_pd
                     (pid_control_V1_B.u4, pid_control_V1_B.u5));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp);
    pid_control_V1_B.u5 = pid_control_V1_B.dv[0];
    pid_control_V1_B.u4 = pid_control_V1_B.dv[1];
    pid_control_V1_B.Va_b_idx_2 += pid_control_V1_B.c_theta[3 *
      pid_control_V1_B.i + 2] * pid_control_V1_B.CD_iw_IGE[pid_control_V1_B.i];
  }

  tmp = _mm_set_pd(pid_control_V1_B.x[3], pid_control_V1_B.x[5]);
  tmp_2 = _mm_sub_pd(_mm_mul_pd(_mm_set_pd(pid_control_V1_B.x[0],
    pid_control_V1_B.x[2]), _mm_loadu_pd(&pid_control_V1_B.x[4])), _mm_mul_pd
                     (_mm_loadu_pd(&pid_control_V1_B.x[1]), tmp));
  _mm_storeu_pd(&pid_control_V1_B.CD_iw_IGE[0], tmp_2);
  pid_control_V1_B.CD_iw_IGE[2] = pid_control_V1_B.x[1] * pid_control_V1_B.x[3]
    - pid_control_V1_B.x[0] * pid_control_V1_B.x[4];
  tmp = _mm_sub_pd(_mm_set_pd(pid_control_V1_B.Mcg_b_idx_1,
    pid_control_V1_B.c_theta_m), _mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (pid_control_V1_B.Vb_w[0], pid_control_V1_B.Vb_w[2]), _mm_loadu_pd
    (&pid_control_V1_B.x[4])), _mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.Vb_w[1]),
    tmp)));
  _mm_storeu_pd(&pid_control_V1_B.Mcg_b[0], tmp);
  pid_control_V1_B.Mcg_b[2] = pid_control_V1_B.Mcg_b_idx_2 -
    (pid_control_V1_B.Vb_w[1] * pid_control_V1_B.x[3] - pid_control_V1_B.Vb_w[0]
     * pid_control_V1_B.x[4]);
  pid_control_V1_B.c_theta_tmp_g1 = 0.0;
  pid_control_V1_B.c_theta_tmp_g = 0.0;
  pid_control_V1_B.c_theta_tmp = 0.0;
  for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 3; pid_control_V1_B.i++) {
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
      (&a[3 * pid_control_V1_B.i]), _mm_set1_pd
      (pid_control_V1_B.Mcg_b[pid_control_V1_B.i])), _mm_set_pd
      (pid_control_V1_B.c_theta_tmp_g, pid_control_V1_B.c_theta_tmp_g1)));
    pid_control_V1_B.c_theta_tmp_g1 = pid_control_V1_B.dv[0];
    pid_control_V1_B.c_theta_tmp_g = pid_control_V1_B.dv[1];
    pid_control_V1_B.c_theta_tmp += a_0[3 * pid_control_V1_B.i + 2] *
      pid_control_V1_B.Mcg_b[pid_control_V1_B.i];
    pid_control_V1_B.c_theta[3 * pid_control_V1_B.i] = c[pid_control_V1_B.i];
  }

  pid_control_V1_B.Vb_w[2] = pid_control_V1_B.c_theta_tmp;
  pid_control_V1_B.Vb_w[1] = pid_control_V1_B.c_theta_tmp_g;
  pid_control_V1_B.Vb_w[0] = pid_control_V1_B.c_theta_tmp_g1;
  pid_control_V1_B.c_theta[1] = 0.0;
  pid_control_V1_B.c_theta[4] = pid_control_V1_B.FA_b_tmp_k;
  pid_control_V1_B.c_theta[7] = -sin(pid_control_V1_B.s_psi);
  pid_control_V1_B.c_theta[2] = 0.0;
  pid_control_V1_B.c_theta[5] = pid_control_V1_B.FA_b_tmp;
  pid_control_V1_B.c_theta[8] = pid_control_V1_B.FA_b_tmp_k;
  pid_control_V1_B.c_theta_tmp_g1 = 0.0;
  pid_control_V1_B.c_theta_tmp_g = 0.0;
  pid_control_V1_B.c_theta_tmp = 0.0;
  for (pid_control_V1_B.i = 0; pid_control_V1_B.i < 3; pid_control_V1_B.i++) {
    tmp = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.c_theta[3 *
      pid_control_V1_B.i]), _mm_set1_pd
      (pid_control_V1_B.wbe_b[pid_control_V1_B.i])), _mm_set_pd
                     (pid_control_V1_B.c_theta_tmp_g,
                      pid_control_V1_B.c_theta_tmp_g1));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp);
    pid_control_V1_B.c_theta_tmp_g1 = pid_control_V1_B.dv[0];
    pid_control_V1_B.c_theta_tmp_g = pid_control_V1_B.dv[1];
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (0.0089285714285714281, pid_control_V1_B.c_theta[3 * pid_control_V1_B.i +
       2]), _mm_set_pd(pid_control_V1_B.F_b[pid_control_V1_B.i],
                       pid_control_V1_B.wbe_b[pid_control_V1_B.i])), _mm_mul_pd
      (_mm_set_pd(pid_control_V1_B.CD_iw_IGE[pid_control_V1_B.i],
                  pid_control_V1_B.c_theta_tmp), _mm_set_pd(-1.0, 1.0))));
    pid_control_V1_B.c_theta_tmp = pid_control_V1_B.dv[0];
    pid_control_V1_B.XDOT[pid_control_V1_B.i] = pid_control_V1_B.dv[1];
    pid_control_V1_B.XDOT[pid_control_V1_B.i + 3] =
      pid_control_V1_B.Vb_w[pid_control_V1_B.i];
  }

  pid_control_V1_B.XDOT[9] = pid_control_V1_B.u5;
  pid_control_V1_B.XDOT[10] = pid_control_V1_B.u4;
  pid_control_V1_B.XDOT[11] = pid_control_V1_B.Va_b_idx_2;
  pid_control_V1_B.XDOT[12] = pid_control_V1_B.CL_w_IGE /
    (pid_control_V1_B.CD_iw_IGE_c + 0.0306);
  pid_control_V1_B.XDOT[19] = 0.0;
  pid_control_V1_B.XDOT[20] = pid_control_V1_B.Cl;
  pid_control_V1_B.XDOT[21] = pid_control_V1_B.Cm;
  pid_control_V1_B.XDOT[22] = pid_control_V1_B.Cn;
  pid_control_V1_B.XDOT[23] = pid_control_V1_B.c_psi;
  pid_control_V1_B.XDOT[24] = pid_control_V1_B.s_psi;
  pid_control_V1_B.XDOT[25] = pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.XDOT[26] = pid_control_V1_B.CL_h_OGE;
  pid_control_V1_B.XDOT[27] = pid_control_V1_B.CL_w_IGE;
  pid_control_V1_B.XDOT[28] = pid_control_V1_B.u2;
  pid_control_V1_B.XDOT[29] = pid_control_V1_B.CD_iw_IGE_c;
  pid_control_V1_B.XDOT[30] = pid_control_V1_B.CD_ih_IGE;
  pid_control_V1_B.XDOT[6] = pid_control_V1_B.c_theta_tmp_g1;
  pid_control_V1_B.XDOT[13] = pid_control_V1_B.F_b[0];
  pid_control_V1_B.XDOT[16] = pid_control_V1_B.c_theta_m;
  pid_control_V1_B.XDOT[31] = pid_control_V1_B.Va;
  pid_control_V1_B.XDOT[34] = pid_control_V1_B.s_theta;
  pid_control_V1_B.XDOT[37] = pid_control_V1_B.FA_b_idx_0;
  pid_control_V1_B.XDOT[7] = pid_control_V1_B.c_theta_tmp_g;
  pid_control_V1_B.XDOT[14] = pid_control_V1_B.F_b[1];
  pid_control_V1_B.XDOT[17] = pid_control_V1_B.Mcg_b_idx_1;
  pid_control_V1_B.XDOT[32] = pid_control_V1_B.s_phi;
  pid_control_V1_B.XDOT[35] = 0.0;
  pid_control_V1_B.XDOT[38] = pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.XDOT[8] = pid_control_V1_B.c_theta_tmp;
  pid_control_V1_B.XDOT[15] = pid_control_V1_B.F_b[2];
  pid_control_V1_B.XDOT[18] = pid_control_V1_B.Mcg_b_idx_2;
  pid_control_V1_B.XDOT[33] = pid_control_V1_B.c_phi;
  pid_control_V1_B.XDOT[36] = 0.17431148549531633;
  pid_control_V1_B.XDOT[39] = pid_control_V1_B.FA_b_idx_2;
  if (pid_control_V1_B.b) {
  }

  /* Gain: '<Root>/Gain' */
  pid_control_V1_B.Gain = -pid_control_V1_B.x[11];
  if (pid_control_V1_B.b) {
  }

  /* BusAssignment: '<Root>/Bus Assignment' */
  memset(&pid_control_V1_B.BusAssignment, 0, sizeof
         (SL_Bus_gazebo_msgs_SetEntityStateRequest));

  /* Gain: '<Root>/Gain2' incorporates:
   *  Gain: '<Root>/Gain1'
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   */
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_div_pd(_mm_set_pd
    (-pid_control_V1_B.x[7], -pid_control_V1_B.x[8]), _mm_set1_pd(2.0)));

  /* MATLABSystem: '<Root>/Coordinate Transformation Conversion' incorporates:
   *  Constant: '<Root>/Constant'
   *  Sum: '<Root>/Sum'
   */
  pid_control_V1_B.c_psi = (pid_control_V1_B.x[6] + 1.5707963267948966) / 2.0;
  pid_control_V1_B.s_psi = sin(pid_control_V1_B.dv[0]);
  pid_control_V1_B.CL_w_OGE = sin(pid_control_V1_B.dv[1]);
  pid_control_V1_B.CL_h_OGE = sin(pid_control_V1_B.c_psi);
  pid_control_V1_B.CL_w_IGE = cos(pid_control_V1_B.dv[0]);
  pid_control_V1_B.u2 = cos(pid_control_V1_B.dv[1]);
  pid_control_V1_B.CD_iw_IGE_c = cos(pid_control_V1_B.c_psi);

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  Gain: '<Root>/Gain3'
   */
  pid_control_V1_B.BusAssignment.state.pose.position.x = pid_control_V1_B.x[9];
  pid_control_V1_B.BusAssignment.state.pose.position.y = -pid_control_V1_B.x[10];
  pid_control_V1_B.BusAssignment.state.pose.position.z = pid_control_V1_B.Gain;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_B.c_psi = pid_control_V1_B.CL_w_IGE * pid_control_V1_B.u2;

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   * */
  pid_control_V1_B.BusAssignment.state.pose.orientation.w =
    pid_control_V1_B.s_psi * pid_control_V1_B.CL_w_OGE *
    pid_control_V1_B.CL_h_OGE + pid_control_V1_B.c_psi *
    pid_control_V1_B.CD_iw_IGE_c;
  pid_control_V1_B.BusAssignment.state.pose.orientation.z =
    pid_control_V1_B.c_psi * pid_control_V1_B.CL_h_OGE -
    pid_control_V1_B.CD_iw_IGE_c * pid_control_V1_B.s_psi *
    pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.BusAssignment.state.pose.orientation.y =
    pid_control_V1_B.CL_w_IGE * pid_control_V1_B.CD_iw_IGE_c *
    pid_control_V1_B.CL_w_OGE + pid_control_V1_B.u2 * pid_control_V1_B.s_psi *
    pid_control_V1_B.CL_h_OGE;
  pid_control_V1_B.BusAssignment.state.pose.orientation.x = pid_control_V1_B.u2 *
    pid_control_V1_B.CD_iw_IGE_c * pid_control_V1_B.s_psi -
    pid_control_V1_B.CL_w_IGE * pid_control_V1_B.CL_w_OGE *
    pid_control_V1_B.CL_h_OGE;
  memcpy(&pid_control_V1_B.BusAssignment.state.name[0],
         &pid_control_V1_B.stringOut_l[0], sizeof(uint8_T) << 7U);
  memcpy(&pid_control_V1_B.BusAssignment.state.reference_frame[0],
         &pid_control_V1_B.stringOut[0], sizeof(uint8_T) << 7U);
  pid_control_V1_B.BusAssignment.state.name_SL_Info.CurrentLength =
    pid_control_V1_B.lengthOut_e;
  pid_control_V1_B.BusAssignment.state.reference_frame_SL_Info.CurrentLength =
    pid_control_V1_B.lengthOut;

  /* Outputs for Atomic SubSystem: '<Root>/Call Service' */
  /* MATLABSystem: '<S2>/ServiceCaller' */
  pid_control_V1_B.serverAvailableOnTime =
    ServCall_pid_control_V1_326.waitForServer(5.0);
  if (pid_control_V1_B.serverAvailableOnTime) {
    ServCall_pid_control_V1_326.call(&pid_control_V1_B.BusAssignment,
      &pid_control_V1_B.r);
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
  pid_control_V1_B.Power = sqrt((pid_control_V1_B.x[0] * pid_control_V1_B.x[0] +
    pid_control_V1_B.x[1] * pid_control_V1_B.x[1]) + pid_control_V1_B.x[2] *
    pid_control_V1_B.x[2]) * pid_control_V1_B.XDOT[34];

  /* Gain: '<S7>/Gain3' */
  pid_control_V1_B.Gain3 = 0.001 * pid_control_V1_B.Power;
  if (pid_control_V1_B.b) {
  }

  /* Product: '<S7>/Divide' incorporates:
   *  Constant: '<S7>/thrust efficiency Cp?'
   */
  pid_control_V1_B.powerdemand = pid_control_V1_B.Gain3 / 0.248;
  if (pid_control_V1_B.b) {
  }

  /* Product: '<S7>/Divide1' */
  pid_control_V1_B.loadtorque = pid_control_V1_B.powerdemand /
    pid_control_V1_ConstB.motorspeed;
  if (pid_control_V1_B.b) {
  }

  /* Gain: '<S88>/ZeroGain' */
  pid_control_V1_B.s_psi = 0.0 * pid_control_V1_B.SignPreSat_o;

  /* DeadZone: '<S90>/DeadZone' */
  if (pid_control_V1_B.SignPreSat_o > 1.0) {
    pid_control_V1_B.SignPreSat_o--;
  } else if (pid_control_V1_B.SignPreSat_o >= 0.0) {
    pid_control_V1_B.SignPreSat_o = 0.0;
  }

  /* End of DeadZone: '<S90>/DeadZone' */

  /* Gain: '<S95>/Integral Gain' incorporates:
   *  Constant: '<Root>/Constant2'
   *  Gain: '<Root>/Gain4'
   *  Sum: '<Root>/Sum2'
   */
  pid_control_V1_B.Switch = (1.0 - (-pid_control_V1_B.x[11])) * 0.05;

  /* Signum: '<S88>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat_o)) {
    /* DataTypeConversion: '<S88>/DataTypeConv1' */
    pid_control_V1_B.i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat_o < 0.0) {
      /* DataTypeConversion: '<S88>/DataTypeConv1' */
      pid_control_V1_B.c_psi = -1.0;
    } else {
      /* DataTypeConversion: '<S88>/DataTypeConv1' */
      pid_control_V1_B.c_psi = (pid_control_V1_B.SignPreSat_o > 0.0);
    }

    /* DataTypeConversion: '<S88>/DataTypeConv1' */
    pid_control_V1_B.i = static_cast<int32_T>(fmod(pid_control_V1_B.c_psi, 256.0));
  }

  /* End of Signum: '<S88>/SignPreSat' */

  /* Signum: '<S88>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch)) {
    /* DataTypeConversion: '<S88>/DataTypeConv2' */
    pid_control_V1_B.i_m = 0;
  } else {
    if (pid_control_V1_B.Switch < 0.0) {
      /* DataTypeConversion: '<S88>/DataTypeConv2' */
      pid_control_V1_B.c_psi = -1.0;
    } else {
      /* DataTypeConversion: '<S88>/DataTypeConv2' */
      pid_control_V1_B.c_psi = (pid_control_V1_B.Switch > 0.0);
    }

    /* DataTypeConversion: '<S88>/DataTypeConv2' */
    pid_control_V1_B.i_m = static_cast<int32_T>(fmod(pid_control_V1_B.c_psi,
      256.0));
  }

  /* End of Signum: '<S88>/SignPreIntegrator' */

  /* DataTypeConversion: '<S88>/DataTypeConv1' */
  if (pid_control_V1_B.i < 0) {
    pid_control_V1_B.i = static_cast<int8_T>(-static_cast<int8_T>
      (static_cast<uint8_T>(-static_cast<real_T>(pid_control_V1_B.i))));
  }

  /* DataTypeConversion: '<S88>/DataTypeConv2' */
  if (pid_control_V1_B.i_m < 0) {
    pid_control_V1_B.i_m = static_cast<int8_T>(-static_cast<int8_T>(static_cast<
      uint8_T>(-static_cast<real_T>(pid_control_V1_B.i_m))));
  }

  /* Logic: '<S88>/AND3' incorporates:
   *  DataTypeConversion: '<S88>/DataTypeConv1'
   *  DataTypeConversion: '<S88>/DataTypeConv2'
   *  RelationalOperator: '<S88>/Equal1'
   *  RelationalOperator: '<S88>/NotEqual'
   */
  pid_control_V1_B.AND3 = ((pid_control_V1_B.s_psi !=
    pid_control_V1_B.SignPreSat_o) && (pid_control_V1_B.i ==
    pid_control_V1_B.i_m));
  if (pid_control_V1_B.b) {
    /* Memory: '<S88>/Memory' */
    pid_control_V1_B.Memory = pid_control_V1_DW.Memory_PreviousInput;
  }

  /* Switch: '<S88>/Switch' */
  if (pid_control_V1_B.Memory) {
    /* Gain: '<S95>/Integral Gain' incorporates:
     *  Constant: '<S88>/Constant1'
     *  Switch: '<S88>/Switch'
     */
    pid_control_V1_B.Switch = 0.0;
  }

  /* End of Switch: '<S88>/Switch' */

  /* Gain: '<S34>/ZeroGain' */
  pid_control_V1_B.SignPreSat_o = 0.0 * pid_control_V1_B.SignPreSat;

  /* DeadZone: '<S36>/DeadZone' */
  if (pid_control_V1_B.SignPreSat > 0.3490658503988659) {
    pid_control_V1_B.SignPreSat -= 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat >= -0.3490658503988659) {
    pid_control_V1_B.SignPreSat = 0.0;
  } else {
    pid_control_V1_B.SignPreSat -= -0.3490658503988659;
  }

  /* End of DeadZone: '<S36>/DeadZone' */

  /* Gain: '<S41>/Integral Gain' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Sum: '<Root>/Sum1'
   */
  pid_control_V1_B.Switch_a = (0.0087266462599716477 - pid_control_V1_B.x[7]) *
    1.05;

  /* Signum: '<S34>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat)) {
    /* DataTypeConversion: '<S34>/DataTypeConv1' */
    pid_control_V1_B.i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat < 0.0) {
      /* DataTypeConversion: '<S34>/DataTypeConv1' */
      pid_control_V1_B.c_psi = -1.0;
    } else {
      /* DataTypeConversion: '<S34>/DataTypeConv1' */
      pid_control_V1_B.c_psi = (pid_control_V1_B.SignPreSat > 0.0);
    }

    /* DataTypeConversion: '<S34>/DataTypeConv1' */
    pid_control_V1_B.i = static_cast<int32_T>(fmod(pid_control_V1_B.c_psi, 256.0));
  }

  /* End of Signum: '<S34>/SignPreSat' */

  /* Signum: '<S34>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch_a)) {
    /* DataTypeConversion: '<S34>/DataTypeConv2' */
    pid_control_V1_B.i_m = 0;
  } else {
    if (pid_control_V1_B.Switch_a < 0.0) {
      /* DataTypeConversion: '<S34>/DataTypeConv2' */
      pid_control_V1_B.c_psi = -1.0;
    } else {
      /* DataTypeConversion: '<S34>/DataTypeConv2' */
      pid_control_V1_B.c_psi = (pid_control_V1_B.Switch_a > 0.0);
    }

    /* DataTypeConversion: '<S34>/DataTypeConv2' */
    pid_control_V1_B.i_m = static_cast<int32_T>(fmod(pid_control_V1_B.c_psi,
      256.0));
  }

  /* End of Signum: '<S34>/SignPreIntegrator' */

  /* DataTypeConversion: '<S34>/DataTypeConv1' */
  if (pid_control_V1_B.i < 0) {
    pid_control_V1_B.i = static_cast<int8_T>(-static_cast<int8_T>
      (static_cast<uint8_T>(-static_cast<real_T>(pid_control_V1_B.i))));
  }

  /* DataTypeConversion: '<S34>/DataTypeConv2' */
  if (pid_control_V1_B.i_m < 0) {
    pid_control_V1_B.i_m = static_cast<int8_T>(-static_cast<int8_T>(static_cast<
      uint8_T>(-static_cast<real_T>(pid_control_V1_B.i_m))));
  }

  /* Logic: '<S34>/AND3' incorporates:
   *  DataTypeConversion: '<S34>/DataTypeConv1'
   *  DataTypeConversion: '<S34>/DataTypeConv2'
   *  RelationalOperator: '<S34>/Equal1'
   *  RelationalOperator: '<S34>/NotEqual'
   */
  pid_control_V1_B.AND3_j = ((pid_control_V1_B.SignPreSat_o !=
    pid_control_V1_B.SignPreSat) && (pid_control_V1_B.i == pid_control_V1_B.i_m));
  if (pid_control_V1_B.b) {
    /* Memory: '<S34>/Memory' */
    pid_control_V1_B.Memory_n = pid_control_V1_DW.Memory_PreviousInput_d;
  }

  /* Switch: '<S34>/Switch' */
  if (pid_control_V1_B.Memory_n) {
    /* Gain: '<S41>/Integral Gain' incorporates:
     *  Constant: '<S34>/Constant1'
     *  Switch: '<S34>/Switch'
     */
    pid_control_V1_B.Switch_a = 0.0;
  }

  /* End of Switch: '<S34>/Switch' */

  /* Gain: '<S7>/Gain1' incorporates:
   *  Integrator: '<S7>/Integrator1'
   */
  pid_control_V1_B.EnergykWh = 2.7777777777777776E-7 *
    pid_control_V1_X.Integrator1_CSTATE;
  if (pid_control_V1_B.b) {
  }

  if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
    if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
      /* Update for Memory: '<S88>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput = pid_control_V1_B.AND3;

      /* Update for Memory: '<S34>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput_d = pid_control_V1_B.AND3_j;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
    rt_ertODEUpdateContinuousStates(&(&pid_control_V1_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++(&pid_control_V1_M)->Timing.clockTick0)) {
      ++(&pid_control_V1_M)->Timing.clockTickH0;
    }

    (&pid_control_V1_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&pid_control_V1_M)
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
      (&pid_control_V1_M)->Timing.clockTick1++;
      if (!(&pid_control_V1_M)->Timing.clockTick1) {
        (&pid_control_V1_M)->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void pid_control_V1::pid_control_V1_derivatives()
{
  XDot_pid_control_V1_T *_rtXdot;
  _rtXdot = ((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs);

  /* Derivatives for Integrator: '<S7>/Integrator' */
  memcpy(&_rtXdot->Integrator_CSTATE[0], &pid_control_V1_B.XDOT[0], 12U * sizeof
         (real_T));

  /* Derivatives for Integrator: '<S44>/Integrator' */
  _rtXdot->Integrator_CSTATE_p = pid_control_V1_B.Switch_a;

  /* Derivatives for Integrator: '<S39>/Filter' */
  _rtXdot->Filter_CSTATE = pid_control_V1_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S98>/Integrator' */
  _rtXdot->Integrator_CSTATE_b = pid_control_V1_B.Switch;

  /* Derivatives for Integrator: '<S93>/Filter' */
  _rtXdot->Filter_CSTATE_f = pid_control_V1_B.FilterCoefficient_g;

  /* Derivatives for Integrator: '<S7>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = pid_control_V1_B.Power;
}

/* Model initialize function */
void pid_control_V1::initialize()
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&pid_control_V1_M)->solverInfo, &(&pid_control_V1_M
                          )->Timing.simTimeStep);
    rtsiSetTPtr(&(&pid_control_V1_M)->solverInfo, &rtmGetTPtr((&pid_control_V1_M)));
    rtsiSetStepSizePtr(&(&pid_control_V1_M)->solverInfo, &(&pid_control_V1_M)
                       ->Timing.stepSize0);
    rtsiSetdXPtr(&(&pid_control_V1_M)->solverInfo, &(&pid_control_V1_M)->derivs);
    rtsiSetContStatesPtr(&(&pid_control_V1_M)->solverInfo, (real_T **)
                         &(&pid_control_V1_M)->contStates);
    rtsiSetNumContStatesPtr(&(&pid_control_V1_M)->solverInfo,
      &(&pid_control_V1_M)->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&pid_control_V1_M)->solverInfo,
      &(&pid_control_V1_M)->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&pid_control_V1_M)->solverInfo,
      &(&pid_control_V1_M)->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&pid_control_V1_M)->solverInfo,
      &(&pid_control_V1_M)->periodicContStateRanges);
    rtsiSetContStateDisabledPtr(&(&pid_control_V1_M)->solverInfo, (boolean_T**)
      &(&pid_control_V1_M)->contStateDisabled);
    rtsiSetErrorStatusPtr(&(&pid_control_V1_M)->solverInfo, (&rtmGetErrorStatus
      ((&pid_control_V1_M))));
    rtsiSetRTModelPtr(&(&pid_control_V1_M)->solverInfo, (&pid_control_V1_M));
  }

  rtsiSetSimTimeStep(&(&pid_control_V1_M)->solverInfo, MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange(&(&pid_control_V1_M)->solverInfo, false);
  rtsiSetIsContModeFrozen(&(&pid_control_V1_M)->solverInfo, false);
  (&pid_control_V1_M)->intgData.y = (&pid_control_V1_M)->odeY;
  (&pid_control_V1_M)->intgData.f[0] = (&pid_control_V1_M)->odeF[0];
  (&pid_control_V1_M)->intgData.f[1] = (&pid_control_V1_M)->odeF[1];
  (&pid_control_V1_M)->intgData.f[2] = (&pid_control_V1_M)->odeF[2];
  (&pid_control_V1_M)->contStates = ((X_pid_control_V1_T *) &pid_control_V1_X);
  (&pid_control_V1_M)->contStateDisabled = ((XDis_pid_control_V1_T *)
    &pid_control_V1_XDis);
  (&pid_control_V1_M)->Timing.tStart = (0.0);
  rtsiSetSolverData(&(&pid_control_V1_M)->solverInfo, static_cast<void *>
                    (&(&pid_control_V1_M)->intgData));
  rtsiSetSolverName(&(&pid_control_V1_M)->solverInfo,"ode3");
  rtmSetTPtr((&pid_control_V1_M), &(&pid_control_V1_M)->Timing.tArray[0]);
  (&pid_control_V1_M)->Timing.stepSize0 = 0.01;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_DW.objisempty = true;
  pid_control_V1_DW.obj_c.isInitialized = 1;

  /* Start for Atomic SubSystem: '<Root>/Call Service' */
  /* Start for MATLABSystem: '<S2>/ServiceCaller' */
  pid_control_V1_DW.obj.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_f = true;
  pid_control_V1_DW.obj.isSetupComplete = false;
  pid_control_V1_DW.obj.isInitialized = 1;
  pid_con_ServiceCaller_setupImpl(&pid_control_V1_DW.obj);
  pid_control_V1_DW.obj.isSetupComplete = true;

  /* End of Start for SubSystem: '<Root>/Call Service' */

  /* InitializeConditions for Integrator: '<S7>/Integrator' */
  memcpy(&pid_control_V1_X.Integrator_CSTATE[0],
         &pid_control_V1_ConstP.Integrator_IC[0], 12U * sizeof(real_T));

  /* InitializeConditions for Integrator: '<S44>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S39>/Filter' */
  pid_control_V1_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S98>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_b = 0.0;

  /* InitializeConditions for Integrator: '<S93>/Filter' */
  pid_control_V1_X.Filter_CSTATE_f = 0.0;

  /* InitializeConditions for Integrator: '<S7>/Integrator1' */
  pid_control_V1_X.Integrator1_CSTATE = 0.0;
}

/* Model terminate function */
void pid_control_V1::terminate()
{
  /* Terminate for Atomic SubSystem: '<Root>/Call Service' */
  /* Terminate for MATLABSystem: '<S2>/ServiceCaller' */
  if (!pid_control_V1_DW.obj.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj.isInitialized == 1) &&
        pid_control_V1_DW.obj.isSetupComplete) {
      ServCall_pid_control_V1_326.resetSvcClientPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S2>/ServiceCaller' */
  /* End of Terminate for SubSystem: '<Root>/Call Service' */
}

/* Constructor */
pid_control_V1::pid_control_V1() :
  pid_control_V1_B(),
  pid_control_V1_DW(),
  pid_control_V1_X(),
  pid_control_V1_XDis(),
  pid_control_V1_M()
{
  /* Currently there is no constructor body generated.*/
}

/* Destructor */
pid_control_V1::~pid_control_V1()
{
  /* Currently there is no destructor body generated.*/
}

/* Real-Time Model get method */
RT_MODEL_pid_control_V1_T * pid_control_V1::getRTM()
{
  return (&pid_control_V1_M);
}
