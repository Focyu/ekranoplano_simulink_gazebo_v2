/*
 * pid_control_V1.cpp
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.37
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Mon Mar  2 19:08:57 2026
 *
 * Target selection: ert.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "pid_control_V1.h"
#include "rtwtypes.h"
#include "pid_control_V1_types.h"
#include <string.h>
#include <math.h>
#include "pid_control_V1_private.h"

extern "C"
{

#include "rt_nonfinite.h"

}

#include <emmintrin.h>
#include "rmw/qos_profiles.h"
#include <stddef.h>
#include "rt_defines.h"

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
  int_T nXc = 21;
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

/*
 * System initialize for enable system:
 *    '<S9>/Enabled Subsystem'
 *    '<S10>/Enabled Subsystem'
 *    '<S11>/Enabled Subsystem'
 */
void pid_control_V1::pid_contr_EnabledSubsystem_Init
  (B_EnabledSubsystem_pid_contro_T *localB)
{
  /* SystemInitialize for SignalConversion generated from: '<S225>/In1' */
  memset(&localB->In1, 0, sizeof(SL_Bus_std_msgs_Float64));
}

/*
 * Output and update for enable system:
 *    '<S9>/Enabled Subsystem'
 *    '<S10>/Enabled Subsystem'
 *    '<S11>/Enabled Subsystem'
 */
void pid_control_V1::pid_control_V1_EnabledSubsystem(boolean_T rtu_Enable, const
  SL_Bus_std_msgs_Float64 *rtu_In1, B_EnabledSubsystem_pid_contro_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S9>/Enabled Subsystem' incorporates:
   *  EnablePort: '<S225>/Enable'
   */
  if (rtu_Enable) {
    /* SignalConversion generated from: '<S225>/In1' */
    localB->In1 = *rtu_In1;
  }

  /* End of Outputs for SubSystem: '<S9>/Enabled Subsystem' */
}

void pid_control_V1::pid_cont_Subscriber_setupImpl_o(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[16] = "/setpoint/pitch";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S10>/SourceBlock' */
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)10.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 16; i++) {
    /* Start for MATLABSystem: '<S10>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_c[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_370.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_c[0],
    qos_profile);
}

void pid_control_V1::pid_con_Subscriber_setupImpl_on(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[14];
  static const char_T b_zeroDelimTopic_0[14] = "/setpoint/yaw";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S11>/SourceBlock' */
  pid_control_V1_B.deadline.sec = 0.0;
  pid_control_V1_B.deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)10.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE, pid_control_V1_B.deadline,
                 lifespan, RMW_QOS_POLICY_LIVELINESS_AUTOMATIC,
                 liveliness_lease_duration, (bool)
                 obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 14; i++) {
    /* Start for MATLABSystem: '<S11>/SourceBlock' */
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_pid_control_V1_377.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void pid_control_V1::pid_contro_Subscriber_setupImpl(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[17] = "/setpoint/altura";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S9>/SourceBlock' */
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)10.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 17; i++) {
    /* Start for MATLABSystem: '<S9>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_m[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_366.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_m[0],
    qos_profile);
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

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int32_T tmp;
    int32_T tmp_0;
    if (u0 > 0.0) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u1 > 0.0) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = atan2(static_cast<real_T>(tmp), static_cast<real_T>(tmp_0));
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

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
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
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
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
  static const real_T a[9] = { 39.71, 0.0, 8.97, 0.0, 85.51, 0.0, 8.97, 0.0,
    114.39 };

  static const real_T b_a[9] = { 0.025636681229753877, -0.0,
    -0.0020103245968257043, -0.0, 0.01169453865045024, -0.0,
    -0.0020103245968257038, -0.0, 0.0088996644080210369 };

  __m128d tmp_2;
  __m128d tmp_3;
  SL_Bus_gazebo_msgs_SetEntityStateResponse tmp;
  int32_T i;
  int32_T tmp_1;
  boolean_T serverAvailableOnTime;
  boolean_T tmp_0;
  static const uint8_T b[11] = { 101U, 107U, 114U, 97U, 110U, 111U, 112U, 108U,
    97U, 110U, 111U };

  static const uint8_T b_0[5] = { 119U, 111U, 114U, 108U, 100U };

  static const real_T a_0[9] = { 39.71, 0.0, 8.97, 0.0, 85.51, 0.0, 8.97, 0.0,
    114.39 };

  static const real_T b_a_0[9] = { 0.025636681229753877, -0.0,
    -0.0020103245968257043, -0.0, 0.01169453865045024, -0.0,
    -0.0020103245968257038, -0.0, 0.0088996644080210369 };

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

  tmp_0 = rtmIsMajorTimeStep((&pid_control_V1_M));
  if (tmp_0) {
    /* MATLAB Function: '<Root>/MATLAB Function' */
    memset(&pid_control_V1_B.stringOut_l[0], 0, sizeof(uint8_T) << 7U);
    for (i = 0; i < 11; i++) {
      pid_control_V1_B.stringOut_l[i] = b[i];
    }

    pid_control_V1_B.lengthOut_e = 11U;

    /* End of MATLAB Function: '<Root>/MATLAB Function' */

    /* MATLAB Function: '<Root>/MATLAB Function1' */
    memset(&pid_control_V1_B.stringOut[0], 0, sizeof(uint8_T) << 7U);
    for (i = 0; i < 5; i++) {
      pid_control_V1_B.stringOut[i] = b_0[i];
    }

    pid_control_V1_B.lengthOut = 5U;

    /* End of MATLAB Function: '<Root>/MATLAB Function1' */
  }

  /* Integrator: '<S12>/Integrator' */
  memcpy(&pid_control_V1_B.x[0], &pid_control_V1_X.Integrator_CSTATE[0], 12U *
         sizeof(real_T));
  if (tmp_0) {
    /* MATLABSystem: '<S10>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_l = Sub_pid_control_V1_370.getLatestMessage(
      &pid_control_V1_B.SourceBlock_o2_e);

    /* Outputs for Enabled SubSystem: '<S10>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_l,
      &pid_control_V1_B.SourceBlock_o2_e, &pid_control_V1_B.EnabledSubsystem_c);

    /* End of Outputs for SubSystem: '<S10>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch1' */
    if (pid_control_V1_B.SourceBlock_o1_l) {
      /* Switch: '<Root>/Switch1' */
      pid_control_V1_B.Switch1 = pid_control_V1_B.EnabledSubsystem_c.In1.data;
    } else {
      /* Switch: '<Root>/Switch1' incorporates:
       *  UnitDelay: '<Root>/Unit Delay1'
       */
      pid_control_V1_B.Switch1 = pid_control_V1_DW.UnitDelay1_DSTATE;
    }

    /* End of Switch: '<Root>/Switch1' */
  }

  /* Gain: '<S98>/Integral Gain' incorporates:
   *  Sum: '<Root>/Sum1'
   */
  pid_control_V1_B.Switch_a = pid_control_V1_B.Switch1 - pid_control_V1_B.x[7];

  /* Gain: '<S104>/Filter Coefficient' incorporates:
   *  Gain: '<S94>/Derivative Gain'
   *  Integrator: '<S96>/Filter'
   *  Sum: '<S96>/SumD'
   */
  pid_control_V1_B.FilterCoefficient = (-12.5 * pid_control_V1_B.Switch_a -
    pid_control_V1_X.Filter_CSTATE) * 100.0;

  /* Sum: '<S110>/Sum' incorporates:
   *  Gain: '<S106>/Proportional Gain'
   *  Integrator: '<S101>/Integrator'
   */
  pid_control_V1_B.SignPreSat = (-25.0 * pid_control_V1_B.Switch_a +
    pid_control_V1_X.Integrator_CSTATE_p) + pid_control_V1_B.FilterCoefficient;

  /* Saturate: '<S108>/Saturation' */
  if (pid_control_V1_B.SignPreSat > 0.3490658503988659) {
    /* Saturate: '<S108>/Saturation' */
    pid_control_V1_B.Saturation = 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat < -0.3490658503988659) {
    /* Saturate: '<S108>/Saturation' */
    pid_control_V1_B.Saturation = -0.3490658503988659;
  } else {
    /* Saturate: '<S108>/Saturation' */
    pid_control_V1_B.Saturation = pid_control_V1_B.SignPreSat;
  }

  /* End of Saturate: '<S108>/Saturation' */
  if (tmp_0) {
    /* MATLABSystem: '<S11>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1 = Sub_pid_control_V1_377.getLatestMessage
      (&pid_control_V1_B.SourceBlock_o2);

    /* Outputs for Enabled SubSystem: '<S11>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1,
      &pid_control_V1_B.SourceBlock_o2, &pid_control_V1_B.EnabledSubsystem_a);

    /* End of Outputs for SubSystem: '<S11>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch2' */
    if (pid_control_V1_B.SourceBlock_o1) {
      /* Switch: '<Root>/Switch2' */
      pid_control_V1_B.Switch2 = pid_control_V1_B.EnabledSubsystem_a.In1.data;
    } else {
      /* Switch: '<Root>/Switch2' incorporates:
       *  UnitDelay: '<Root>/Unit Delay2'
       */
      pid_control_V1_B.Switch2 = pid_control_V1_DW.UnitDelay2_DSTATE;
    }

    /* End of Switch: '<Root>/Switch2' */
  }

  /* Sum: '<Root>/Sum5' */
  pid_control_V1_B.Sum5 = pid_control_V1_B.Switch2 - pid_control_V1_B.x[8];

  /* Gain: '<S210>/Filter Coefficient' incorporates:
   *  Gain: '<S200>/Derivative Gain'
   *  Integrator: '<S202>/Filter'
   *  Sum: '<S202>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_p = (-5.0 * pid_control_V1_B.Sum5 -
    pid_control_V1_X.Filter_CSTATE_f) * 100.0;

  /* Sum: '<S216>/Sum' incorporates:
   *  Gain: '<S212>/Proportional Gain'
   *  Integrator: '<S207>/Integrator'
   */
  pid_control_V1_B.Saturation_m = (-15.0 * pid_control_V1_B.Sum5 +
    pid_control_V1_X.Integrator_CSTATE_d) + pid_control_V1_B.FilterCoefficient_p;

  /* Saturate: '<S214>/Saturation' */
  if (pid_control_V1_B.Saturation_m > 0.26179938779914941) {
    /* Sum: '<S216>/Sum' incorporates:
     *  Saturate: '<S214>/Saturation'
     */
    pid_control_V1_B.Saturation_m = 0.26179938779914941;
  } else if (pid_control_V1_B.Saturation_m < -0.26179938779914941) {
    /* Sum: '<S216>/Sum' incorporates:
     *  Saturate: '<S214>/Saturation'
     */
    pid_control_V1_B.Saturation_m = -0.26179938779914941;
  }

  /* End of Saturate: '<S214>/Saturation' */

  /* Gain: '<S50>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant4'
   *  Gain: '<S40>/Derivative Gain'
   *  Integrator: '<S42>/Filter'
   *  Sum: '<Root>/Sum4'
   *  Sum: '<S42>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_c = ((0.0 - pid_control_V1_B.x[6]) * -10.0
    - pid_control_V1_X.Filter_CSTATE_g) * 100.0;

  /* Sum: '<S56>/Sum' incorporates:
   *  Constant: '<Root>/Constant4'
   *  Gain: '<S52>/Proportional Gain'
   *  Integrator: '<S47>/Integrator'
   *  Sum: '<Root>/Sum4'
   */
  pid_control_V1_B.Saturation_k = ((0.0 - pid_control_V1_B.x[6]) * -20.0 +
    pid_control_V1_X.Integrator_CSTATE_m) + pid_control_V1_B.FilterCoefficient_c;

  /* Saturate: '<S54>/Saturation' */
  if (pid_control_V1_B.Saturation_k > 0.3490658503988659) {
    /* Sum: '<S56>/Sum' incorporates:
     *  Saturate: '<S54>/Saturation'
     */
    pid_control_V1_B.Saturation_k = 0.3490658503988659;
  } else if (pid_control_V1_B.Saturation_k < -0.3490658503988659) {
    /* Sum: '<S56>/Sum' incorporates:
     *  Saturate: '<S54>/Saturation'
     */
    pid_control_V1_B.Saturation_k = -0.3490658503988659;
  }

  /* End of Saturate: '<S54>/Saturation' */
  if (tmp_0) {
    /* MATLABSystem: '<S9>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_i = Sub_pid_control_V1_366.getLatestMessage(
      &pid_control_V1_B.SourceBlock_o2_l);

    /* Outputs for Enabled SubSystem: '<S9>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_i,
      &pid_control_V1_B.SourceBlock_o2_l, &pid_control_V1_B.EnabledSubsystem);

    /* End of Outputs for SubSystem: '<S9>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch' */
    if (pid_control_V1_B.SourceBlock_o1_i) {
      /* Switch: '<Root>/Switch' */
      pid_control_V1_B.Switch = pid_control_V1_B.EnabledSubsystem.In1.data;
    } else {
      /* Switch: '<Root>/Switch' incorporates:
       *  UnitDelay: '<Root>/Unit Delay'
       */
      pid_control_V1_B.Switch = pid_control_V1_DW.UnitDelay_DSTATE;
    }

    /* End of Switch: '<Root>/Switch' */
  }

  /* Gain: '<S152>/Integral Gain' incorporates:
   *  Gain: '<Root>/Gain4'
   *  Sum: '<Root>/Sum2'
   */
  pid_control_V1_B.Switch_c = pid_control_V1_B.Switch - (-pid_control_V1_B.x[11]);

  /* Gain: '<S148>/Derivative Gain' incorporates:
   *  Gain: '<S160>/Proportional Gain'
   */
  pid_control_V1_B.SignPreSat_o = 0.25 * pid_control_V1_B.Switch_c;

  /* Gain: '<S158>/Filter Coefficient' incorporates:
   *  Gain: '<S148>/Derivative Gain'
   *  Integrator: '<S150>/Filter'
   *  Sum: '<S150>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_g = (pid_control_V1_B.SignPreSat_o -
    pid_control_V1_X.Filter_CSTATE_fy) * 100.0;

  /* Sum: '<S164>/Sum' incorporates:
   *  Integrator: '<S155>/Integrator'
   */
  pid_control_V1_B.SignPreSat_o = (pid_control_V1_B.SignPreSat_o +
    pid_control_V1_X.Integrator_CSTATE_b) + pid_control_V1_B.FilterCoefficient_g;

  /* Saturate: '<S162>/Saturation' */
  if (pid_control_V1_B.SignPreSat_o > 1.0) {
    /* Saturate: '<S162>/Saturation' */
    pid_control_V1_B.Saturation_p = 1.0;
  } else if (pid_control_V1_B.SignPreSat_o < 0.0) {
    /* Saturate: '<S162>/Saturation' */
    pid_control_V1_B.Saturation_p = 0.0;
  } else {
    /* Saturate: '<S162>/Saturation' */
    pid_control_V1_B.Saturation_p = pid_control_V1_B.SignPreSat_o;
  }

  /* End of Saturate: '<S162>/Saturation' */

  /* Saturate: '<Root>/Saturation' incorporates:
   *  Constant: '<Root>/Constant3'
   *  Sum: '<Root>/Sum3'
   */
  if (pid_control_V1_B.Saturation_p + 0.35 > 1.0) {
    pid_control_V1_B.Vd1 = 1.0;
  } else if (pid_control_V1_B.Saturation_p + 0.35 < 0.3) {
    pid_control_V1_B.Vd1 = 0.3;
  } else {
    pid_control_V1_B.Vd1 = pid_control_V1_B.Saturation_p + 0.35;
  }

  /* Gain: '<Root>/Gain5' incorporates:
   *  Saturate: '<Root>/Saturation'
   */
  pid_control_V1_B.Gain5 = 0.000125 * pid_control_V1_B.Vd1;

  /* SignalConversion generated from: '<S228>/ SFunction ' incorporates:
   *  MATLAB Function: '<S12>/MATLAB Function'
   */
  pid_control_V1_B.TmpSignalConversionAtSFunct[0] =
    pid_control_V1_B.Saturation_k;
  pid_control_V1_B.TmpSignalConversionAtSFunct[1] = pid_control_V1_B.Saturation;
  pid_control_V1_B.TmpSignalConversionAtSFunct[2] =
    pid_control_V1_B.Saturation_m;
  pid_control_V1_B.TmpSignalConversionAtSFunct[3] = pid_control_V1_B.Gain5;
  pid_control_V1_B.TmpSignalConversionAtSFunct[4] = pid_control_V1_B.Gain5;

  /* MATLAB Function: '<S12>/MATLAB Function' */
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[1] <= 0.3490658503988659) {
    pid_control_V1_B.u2 = pid_control_V1_B.TmpSignalConversionAtSFunct[1];
  } else {
    pid_control_V1_B.u2 = 0.3490658503988659;
  }

  if (!(pid_control_V1_B.u2 >= -0.3490658503988659)) {
    pid_control_V1_B.u2 = -0.3490658503988659;
  }

  pid_control_V1_B.Va = sqrt((pid_control_V1_B.x[0] * pid_control_V1_B.x[0] +
    pid_control_V1_B.x[1] * pid_control_V1_B.x[1]) + pid_control_V1_B.x[2] *
    pid_control_V1_B.x[2]);
  if (pid_control_V1_B.Va == 0.0) {
    pid_control_V1_B.Va = 0.001;
  }

  pid_control_V1_B.alpha = rt_atan2d_snf(pid_control_V1_B.x[2],
    pid_control_V1_B.x[0]);
  pid_control_V1_B.beta = asin(pid_control_V1_B.x[1] / pid_control_V1_B.Va);
  if ((-pid_control_V1_B.x[11] - 0.363 <= 0.001) || rtIsNaN(-pid_control_V1_B.x
       [11] - 0.363)) {
    pid_control_V1_B.hw = 0.001;
  } else {
    pid_control_V1_B.hw = -pid_control_V1_B.x[11] - 0.363;
  }

  pid_control_V1_B.Q_k = pid_control_V1_B.Va * pid_control_V1_B.Va * 0.6125;
  pid_control_V1_B.wbe_b[0] = pid_control_V1_B.x[3];
  pid_control_V1_B.wbe_b[1] = pid_control_V1_B.x[4];
  pid_control_V1_B.wbe_b[2] = pid_control_V1_B.x[5];
  pid_control_V1_B.CL_w_OGE = ((pid_control_V1_B.alpha - -0.065449846949787352)
    + 0.043633231299858237) * 4.9604094530365153;
  pid_control_V1_B.CL_h_OGE = ((pid_control_V1_B.alpha - -0.074176493209759012)
    + 0.026179938779914941) * 4.8387748917360032;
  pid_control_V1_B.Cn = pid_control_V1_B.hw / 5.02;
  pid_control_V1_B.hw = (rt_powd_snf(pid_control_V1_B.Cn, 0.787) * 288.0 * exp
    (rt_powd_snf(pid_control_V1_B.Cn, 0.327) * -9.14) * 0.97986308862072491 /
    5.9129476540958859 + 1.0) * pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.FA_b_idx_0 = pid_control_V1_B.u2 * pid_control_V1_B.u2;
  pid_control_V1_B.CL_h_IGE = (rt_powd_snf(fabs((-pid_control_V1_B.x[11] + 0.72)
    / 2.74), 0.787) * 288.0 * exp(rt_powd_snf(fabs((-pid_control_V1_B.x[11] +
    0.72) / 2.74), 0.327) * -9.14) * 0.95628590200128227 / 5.35300902982722 +
    1.0) * ((pid_control_V1_B.FA_b_idx_0 * -0.00141 + 0.0307 *
             pid_control_V1_B.u2) + pid_control_V1_B.CL_h_OGE);
  pid_control_V1_B.CD_iw_IGE = (1.0 - exp(rt_powd_snf(pid_control_V1_B.Cn, 0.686)
    * -10.1)) * (pid_control_V1_B.CL_w_OGE * pid_control_V1_B.CL_w_OGE /
                 21.205750411731103);
  pid_control_V1_B.CD_ih_IGE = (1.0 - exp(rt_powd_snf(fabs((-pid_control_V1_B.x
    [11] + 0.72) / 2.74), 0.686) * -10.1)) * (pid_control_V1_B.CL_h_OGE *
    pid_control_V1_B.CL_h_OGE / 18.943803701146454);
  pid_control_V1_B.CQ = -0.019 * pid_control_V1_B.beta * 180.0 /
    3.1415926535897931;
  pid_control_V1_B.FA_b_idx_1 = sin(pid_control_V1_B.alpha);
  pid_control_V1_B.FA_b_idx_2 = cos(pid_control_V1_B.alpha);
  pid_control_V1_B.FA_b_tmp[0] = pid_control_V1_B.FA_b_idx_2;
  pid_control_V1_B.FA_b_tmp[3] = 0.0;
  pid_control_V1_B.FA_b_tmp[6] = -pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.FA_b_tmp[2] = pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.FA_b_tmp[5] = 0.0;
  pid_control_V1_B.FA_b_tmp[8] = pid_control_V1_B.FA_b_idx_2;
  pid_control_V1_B.Q[0] = -(((pid_control_V1_B.FA_b_idx_0 * -1.08E-5 + 0.000715 *
    pid_control_V1_B.u2) * 1.128 + ((pid_control_V1_B.CD_iw_IGE * 3.334 +
    0.1020204) + pid_control_V1_B.CD_ih_IGE * 1.128)) * pid_control_V1_B.Q_k);
  pid_control_V1_B.Q[1] = pid_control_V1_B.CQ * pid_control_V1_B.Q_k * 3.334;
  pid_control_V1_B.Q[2] = -((pid_control_V1_B.hw * 3.334 +
    pid_control_V1_B.CL_h_IGE * 1.128) * pid_control_V1_B.Q_k);
  pid_control_V1_B.FA_b_tmp[1] = 0.0;
  pid_control_V1_B.FA_b_idx_0 = 0.0;
  pid_control_V1_B.FA_b_tmp[4] = 1.0;
  pid_control_V1_B.FA_b_idx_1 = 0.0;
  pid_control_V1_B.FA_b_tmp[7] = 0.0;
  pid_control_V1_B.FA_b_idx_2 = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_3 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.Q[i])), _mm_set_pd
                       (pid_control_V1_B.FA_b_idx_1, pid_control_V1_B.FA_b_idx_0));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_3);
    pid_control_V1_B.FA_b_idx_0 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FA_b_idx_1 = pid_control_V1_B.dv[1];
    pid_control_V1_B.FA_b_idx_2 += pid_control_V1_B.FA_b_tmp[3 * i + 2] *
      pid_control_V1_B.Q[i];
  }

  if (pid_control_V1_B.TmpSignalConversionAtSFunct[0] <= 0.26179938779914941) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[0];
  } else {
    pid_control_V1_B.Vd1 = 0.26179938779914941;
  }

  if (!(pid_control_V1_B.Vd1 >= -0.3490658503988659)) {
    pid_control_V1_B.Vd1 = -0.3490658503988659;
  }

  pid_control_V1_B.Tp1 = 2.0 * pid_control_V1_B.Va;
  pid_control_V1_B.Cl = (-0.0005 * pid_control_V1_B.beta * 180.0 /
    3.1415926535897931 + pid_control_V1_B.x[3] * 5.02 / pid_control_V1_B.Tp1 *
    -2.0) + -0.5 * pid_control_V1_B.Vd1;
  pid_control_V1_B.u2 = ((exp(pid_control_V1_B.Cn * -4.0) * -0.05 + -1.14 *
    pid_control_V1_B.alpha) + pid_control_V1_B.x[4] * 0.646 /
    pid_control_V1_B.Tp1 * -5.0) + -3.0 * pid_control_V1_B.u2;
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[2] <= 0.26179938779914941) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[2];
  } else {
    pid_control_V1_B.Vd1 = 0.26179938779914941;
  }

  if (!(pid_control_V1_B.Vd1 >= -0.26179938779914941)) {
    pid_control_V1_B.Vd1 = -0.26179938779914941;
  }

  pid_control_V1_B.Cn = (-0.002 * pid_control_V1_B.beta * 180.0 /
    3.1415926535897931 + pid_control_V1_B.x[5] * 5.02 / pid_control_V1_B.Tp1 *
    -1.5) + -0.3 * pid_control_V1_B.Vd1;
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[3] <= 0.000125) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[3];
  } else {
    pid_control_V1_B.Vd1 = 0.000125;
  }

  pid_control_V1_B.Vd1 = pid_control_V1_B.Vd1 / 0.000125 * (37.42 -
    pid_control_V1_B.Va) + pid_control_V1_B.Va;
  pid_control_V1_B.Tp1 = 0.11056515539994617 * pid_control_V1_B.Vd1 *
    (pid_control_V1_B.Vd1 - pid_control_V1_B.Va);
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[4] <= 0.000125) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[4];
  } else {
    pid_control_V1_B.Vd1 = 0.000125;
  }

  pid_control_V1_B.Vd1 = pid_control_V1_B.Vd1 / 0.000125 * (37.42 -
    pid_control_V1_B.Va) + pid_control_V1_B.Va;
  pid_control_V1_B.Vd1 = 0.11056515539994617 * pid_control_V1_B.Vd1 *
    (pid_control_V1_B.Vd1 - pid_control_V1_B.Va);
  pid_control_V1_B.FE1_b_idx_0 = pid_control_V1_B.Tp1 * 0.99619469809174555;
  pid_control_V1_B.FE1_b_idx_1 = 0.0;
  pid_control_V1_B.FE1_b_idx_2 = pid_control_V1_B.Tp1 * 0.087155742747658166;
  pid_control_V1_B.FE2_b[0] = pid_control_V1_B.Vd1 * 0.99619469809174555;
  pid_control_V1_B.FE2_b[2] = pid_control_V1_B.Vd1 * 0.087155742747658166;
  pid_control_V1_B.Vd1 = -9.81 * sin(pid_control_V1_B.x[7]) * 112.0;
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_mul_pd(_mm_mul_pd(_mm_mul_pd
    (_mm_set1_pd(9.81), _mm_set1_pd(cos(pid_control_V1_B.x[7]))), _mm_set_pd(cos
    (pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6]))), _mm_set1_pd(112.0)));
  pid_control_V1_B.Va = pid_control_V1_B.dv[0];
  pid_control_V1_B.Tp1 = pid_control_V1_B.dv[1];
  pid_control_V1_B.Mcg_b_idx_2 = 5.02 * pid_control_V1_B.Q_k * 3.334;
  pid_control_V1_B.Mcg_b_idx_0 = (0.6 * pid_control_V1_B.FE1_b_idx_2 + -0.6 *
    pid_control_V1_B.FE2_b[2]) + pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cl;
  pid_control_V1_B.Q_k = 0.646 * pid_control_V1_B.Q_k * 3.334 *
    pid_control_V1_B.u2 + ((-0.285 * pid_control_V1_B.FE1_b_idx_0 - 0.519 *
    pid_control_V1_B.FE1_b_idx_2) + (-0.285 * pid_control_V1_B.FE2_b[0] - 0.519 *
    pid_control_V1_B.FE2_b[2]));
  pid_control_V1_B.Mcg_b_idx_2 = ((0.0 - 0.6 * pid_control_V1_B.FE1_b_idx_0) +
    (0.0 - -0.6 * pid_control_V1_B.FE2_b[0])) + pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cn;
  pid_control_V1_B.FE_b = pid_control_V1_B.FE1_b_idx_0 + pid_control_V1_B.FE2_b
    [0];
  pid_control_V1_B.FE_b_idx_0 = pid_control_V1_B.FE_b;
  pid_control_V1_B.F_b[0] = (pid_control_V1_B.Vd1 + pid_control_V1_B.FE_b) +
    pid_control_V1_B.FA_b_idx_0;
  pid_control_V1_B.FE1_b_idx_0 = 0.0;
  pid_control_V1_B.F_b[1] = pid_control_V1_B.dv[0] + pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.FE_b = pid_control_V1_B.FE1_b_idx_2 + pid_control_V1_B.FE2_b
    [2];
  pid_control_V1_B.F_b[2] = (pid_control_V1_B.dv[1] + pid_control_V1_B.FE_b) +
    pid_control_V1_B.FA_b_idx_2;
  pid_control_V1_B.FE1_b_idx_2 = 0.0;
  pid_control_V1_B.FA_b_tmp[0] = cos(pid_control_V1_B.x[7]) * cos
    (pid_control_V1_B.x[8]);
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd
    (_mm_set_pd(cos(pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6])),
     _mm_set1_pd(sin(pid_control_V1_B.x[7]))), _mm_set1_pd(cos
    (pid_control_V1_B.x[8]))), _mm_mul_pd(_mm_mul_pd(_mm_set_pd(sin
    (pid_control_V1_B.x[6]), cos(pid_control_V1_B.x[6])), _mm_set1_pd(sin
    (pid_control_V1_B.x[8]))), _mm_set_pd(1.0, -1.0))));
  pid_control_V1_B.FA_b_tmp[3] = pid_control_V1_B.dv[0];
  pid_control_V1_B.FA_b_tmp[6] = pid_control_V1_B.dv[1];
  pid_control_V1_B.FA_b_tmp[1] = cos(pid_control_V1_B.x[7]) * sin
    (pid_control_V1_B.x[8]);
  tmp_3 = _mm_set_pd(-1.0, 1.0);

  /* MATLAB Function: '<S12>/MATLAB Function' */
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd
    (_mm_set_pd(cos(pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6])),
     _mm_set1_pd(sin(pid_control_V1_B.x[7]))), _mm_set1_pd(sin
    (pid_control_V1_B.x[8]))), _mm_mul_pd(_mm_mul_pd(_mm_set_pd(sin
    (pid_control_V1_B.x[6]), cos(pid_control_V1_B.x[6])), _mm_set1_pd(cos
    (pid_control_V1_B.x[8]))), tmp_3)));
  pid_control_V1_B.FA_b_tmp[4] = pid_control_V1_B.dv[0];
  pid_control_V1_B.FA_b_tmp[7] = pid_control_V1_B.dv[1];
  pid_control_V1_B.FA_b_tmp[2] = -sin(pid_control_V1_B.x[7]);
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_mul_pd(_mm_set_pd(cos
    (pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6])), _mm_set1_pd(cos
    (pid_control_V1_B.x[7]))));
  pid_control_V1_B.FA_b_tmp[5] = pid_control_V1_B.dv[0];
  pid_control_V1_B.FA_b_tmp[8] = pid_control_V1_B.dv[1];
  pid_control_V1_B.Q[0] = pid_control_V1_B.x[0];
  pid_control_V1_B.Q[1] = pid_control_V1_B.x[1];
  pid_control_V1_B.Q[2] = pid_control_V1_B.x[2];
  for (i = 0; i < 3; i++) {
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
      (&a[3 * i]), _mm_set1_pd(pid_control_V1_B.wbe_b[i])), _mm_set_pd
      (pid_control_V1_B.FE1_b_idx_1, pid_control_V1_B.FE1_b_idx_0)));
    pid_control_V1_B.FE1_b_idx_0 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FE1_b_idx_1 = pid_control_V1_B.dv[1];
    pid_control_V1_B.FE1_b_idx_2 += a_0[3 * i + 2] * pid_control_V1_B.wbe_b[i];
  }

  pid_control_V1_B.FE2_b_c = 0.0;
  pid_control_V1_B.FE2_b_b = 0.0;
  pid_control_V1_B.FE2_b_p = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.Q[i])), _mm_set_pd(pid_control_V1_B.FE2_b_b,
      pid_control_V1_B.FE2_b_c));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
    pid_control_V1_B.FE2_b_c = pid_control_V1_B.dv[0];
    pid_control_V1_B.FE2_b_b = pid_control_V1_B.dv[1];
    pid_control_V1_B.FE2_b_p += pid_control_V1_B.FA_b_tmp[3 * i + 2] *
      pid_control_V1_B.Q[i];
  }

  tmp_2 = _mm_sub_pd(_mm_mul_pd(_mm_set_pd(pid_control_V1_B.x[0],
    pid_control_V1_B.x[2]), _mm_loadu_pd(&pid_control_V1_B.x[4])), _mm_mul_pd
                     (_mm_loadu_pd(&pid_control_V1_B.x[1]), _mm_set_pd
                      (pid_control_V1_B.x[3], pid_control_V1_B.x[5])));
  _mm_storeu_pd(&pid_control_V1_B.Q[0], tmp_2);
  pid_control_V1_B.Q[2] = pid_control_V1_B.x[1] * pid_control_V1_B.x[3] -
    pid_control_V1_B.x[0] * pid_control_V1_B.x[4];
  _mm_storeu_pd(&pid_control_V1_B.FE2_b[0], _mm_sub_pd(_mm_set_pd
    (pid_control_V1_B.Q_k, pid_control_V1_B.Mcg_b_idx_0), _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(pid_control_V1_B.FE1_b_idx_0, pid_control_V1_B.x[4]), _mm_set_pd
     (pid_control_V1_B.x[5], pid_control_V1_B.FE1_b_idx_2)), _mm_mul_pd
    (_mm_set_pd(pid_control_V1_B.x[3], pid_control_V1_B.FE1_b_idx_1), _mm_set_pd
     (pid_control_V1_B.FE1_b_idx_2, pid_control_V1_B.x[5])))));
  pid_control_V1_B.FE2_b[2] = pid_control_V1_B.Mcg_b_idx_2 -
    (pid_control_V1_B.x[3] * pid_control_V1_B.FE1_b_idx_1 -
     pid_control_V1_B.FE1_b_idx_0 * pid_control_V1_B.x[4]);
  pid_control_V1_B.FE1_b_idx_0 = 0.0;
  pid_control_V1_B.FE1_b_idx_1 = 0.0;
  pid_control_V1_B.FE1_b_idx_2 = 0.0;
  for (i = 0; i < 3; i++) {
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
      (&b_a[3 * i]), _mm_set1_pd(pid_control_V1_B.FE2_b[i])), _mm_set_pd
      (pid_control_V1_B.FE1_b_idx_1, pid_control_V1_B.FE1_b_idx_0)));
    pid_control_V1_B.FE1_b_idx_0 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FE1_b_idx_1 = pid_control_V1_B.dv[1];
    pid_control_V1_B.FE1_b_idx_2 += b_a_0[3 * i + 2] * pid_control_V1_B.FE2_b[i];
  }

  pid_control_V1_B.FE2_b[2] = pid_control_V1_B.FE1_b_idx_2;
  pid_control_V1_B.FE2_b[1] = pid_control_V1_B.FE1_b_idx_1;
  pid_control_V1_B.FE2_b[0] = pid_control_V1_B.FE1_b_idx_0;
  pid_control_V1_B.FA_b_tmp[0] = 1.0;
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_mul_pd(_mm_set_pd(cos
    (pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6])), _mm_set1_pd(tan
    (pid_control_V1_B.x[7]))));
  pid_control_V1_B.FA_b_tmp[3] = pid_control_V1_B.dv[0];
  pid_control_V1_B.FA_b_tmp[6] = pid_control_V1_B.dv[1];
  pid_control_V1_B.FA_b_tmp[1] = 0.0;
  pid_control_V1_B.FA_b_tmp[4] = cos(pid_control_V1_B.x[6]);
  pid_control_V1_B.FA_b_tmp[7] = -sin(pid_control_V1_B.x[6]);
  pid_control_V1_B.FA_b_tmp[2] = 0.0;
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_div_pd(_mm_set_pd(cos
    (pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6])), _mm_set1_pd(cos
    (pid_control_V1_B.x[7]))));
  pid_control_V1_B.FA_b_tmp[5] = pid_control_V1_B.dv[0];
  pid_control_V1_B.FA_b_tmp[8] = pid_control_V1_B.dv[1];
  pid_control_V1_B.FE1_b_idx_0 = 0.0;
  pid_control_V1_B.FE1_b_idx_1 = 0.0;
  pid_control_V1_B.FE1_b_idx_2 = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.wbe_b[i])), _mm_set_pd
                       (pid_control_V1_B.FE1_b_idx_1,
                        pid_control_V1_B.FE1_b_idx_0));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
    pid_control_V1_B.FE1_b_idx_0 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FE1_b_idx_1 = pid_control_V1_B.dv[1];
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (0.0089285714285714281, pid_control_V1_B.FA_b_tmp[3 * i + 2]), _mm_set_pd
      (pid_control_V1_B.F_b[i], pid_control_V1_B.wbe_b[i])), _mm_mul_pd
      (_mm_set_pd(pid_control_V1_B.Q[i], pid_control_V1_B.FE1_b_idx_2), tmp_3)));
    pid_control_V1_B.FE1_b_idx_2 = pid_control_V1_B.dv[0];
    pid_control_V1_B.XDOT[i] = pid_control_V1_B.dv[1];
    pid_control_V1_B.XDOT[i + 3] = pid_control_V1_B.FE2_b[i];
  }

  pid_control_V1_B.XDOT[9] = pid_control_V1_B.FE2_b_c;
  pid_control_V1_B.XDOT[10] = pid_control_V1_B.FE2_b_b;
  pid_control_V1_B.XDOT[11] = pid_control_V1_B.FE2_b_p;
  pid_control_V1_B.XDOT[12] = 0.0;
  pid_control_V1_B.XDOT[19] = pid_control_V1_B.CQ;
  pid_control_V1_B.XDOT[20] = pid_control_V1_B.Cl;
  pid_control_V1_B.XDOT[21] = pid_control_V1_B.u2;
  pid_control_V1_B.XDOT[22] = pid_control_V1_B.Cn;
  pid_control_V1_B.XDOT[23] = pid_control_V1_B.alpha;
  pid_control_V1_B.XDOT[24] = pid_control_V1_B.beta;
  pid_control_V1_B.XDOT[25] = pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.XDOT[26] = pid_control_V1_B.CL_h_OGE;
  pid_control_V1_B.XDOT[27] = pid_control_V1_B.hw;
  pid_control_V1_B.XDOT[28] = pid_control_V1_B.CL_h_IGE;
  pid_control_V1_B.XDOT[29] = pid_control_V1_B.CD_iw_IGE;
  pid_control_V1_B.XDOT[30] = pid_control_V1_B.CD_ih_IGE;
  pid_control_V1_B.XDOT[6] = pid_control_V1_B.FE1_b_idx_0;
  pid_control_V1_B.XDOT[13] = pid_control_V1_B.F_b[0];
  pid_control_V1_B.XDOT[16] = pid_control_V1_B.Mcg_b_idx_0;
  pid_control_V1_B.XDOT[31] = pid_control_V1_B.Vd1;
  pid_control_V1_B.XDOT[34] = pid_control_V1_B.FE_b_idx_0;
  pid_control_V1_B.XDOT[37] = pid_control_V1_B.FA_b_idx_0;
  pid_control_V1_B.XDOT[7] = pid_control_V1_B.FE1_b_idx_1;
  pid_control_V1_B.XDOT[14] = pid_control_V1_B.F_b[1];
  pid_control_V1_B.XDOT[17] = pid_control_V1_B.Q_k;
  pid_control_V1_B.XDOT[32] = pid_control_V1_B.Va;
  pid_control_V1_B.XDOT[35] = 0.0;
  pid_control_V1_B.XDOT[38] = pid_control_V1_B.FA_b_idx_1;
  pid_control_V1_B.XDOT[8] = pid_control_V1_B.FE1_b_idx_2;
  pid_control_V1_B.XDOT[15] = pid_control_V1_B.F_b[2];
  pid_control_V1_B.XDOT[18] = pid_control_V1_B.Mcg_b_idx_2;
  pid_control_V1_B.XDOT[33] = pid_control_V1_B.Tp1;
  pid_control_V1_B.XDOT[36] = pid_control_V1_B.FE_b;
  pid_control_V1_B.XDOT[39] = pid_control_V1_B.FA_b_idx_2;
  if (tmp_0) {
  }

  /* Gain: '<Root>/Gain' */
  pid_control_V1_B.Gain = -pid_control_V1_B.x[11];
  if (tmp_0) {
    /* SignalConversion generated from: '<S12>/To Workspace1' */
    pid_control_V1_B.TmpSignalConversionAtSFunct[0] =
      pid_control_V1_B.Saturation_k;
    pid_control_V1_B.TmpSignalConversionAtSFunct[1] =
      pid_control_V1_B.Saturation;
    pid_control_V1_B.TmpSignalConversionAtSFunct[2] =
      pid_control_V1_B.Saturation_m;
    pid_control_V1_B.TmpSignalConversionAtSFunct[3] = pid_control_V1_B.Gain5;
    pid_control_V1_B.TmpSignalConversionAtSFunct[4] = pid_control_V1_B.Gain5;
  }

  /* BusAssignment: '<Root>/Bus Assignment' */
  memset(&pid_control_V1_B.BusAssignment, 0, sizeof
         (SL_Bus_gazebo_msgs_SetEntityStateRequest));

  /* Gain: '<Root>/Gain2' incorporates:
   *  Gain: '<Root>/Gain1'
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   */
  _mm_storeu_pd(&pid_control_V1_B.wbe_b[0], _mm_div_pd(_mm_set_pd
    (-pid_control_V1_B.x[7], -pid_control_V1_B.x[8]), _mm_set1_pd(2.0)));

  /* MATLABSystem: '<Root>/Coordinate Transformation Conversion' incorporates:
   *  Constant: '<Root>/Constant'
   *  Sum: '<Root>/Sum'
   */
  pid_control_V1_B.wbe_b[2] = (pid_control_V1_B.x[6] + 1.5707963267948966) / 2.0;
  pid_control_V1_B.alpha = sin(pid_control_V1_B.wbe_b[0]);
  pid_control_V1_B.beta = sin(pid_control_V1_B.wbe_b[1]);
  pid_control_V1_B.CL_w_OGE = sin(pid_control_V1_B.wbe_b[2]);
  pid_control_V1_B.CL_h_OGE = cos(pid_control_V1_B.wbe_b[0]);
  pid_control_V1_B.hw = cos(pid_control_V1_B.wbe_b[1]);
  pid_control_V1_B.CL_h_IGE = cos(pid_control_V1_B.wbe_b[2]);

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  Gain: '<Root>/Gain3'
   */
  pid_control_V1_B.BusAssignment.state.pose.position.x = pid_control_V1_B.x[9];
  pid_control_V1_B.BusAssignment.state.pose.position.y = -pid_control_V1_B.x[10];
  pid_control_V1_B.BusAssignment.state.pose.position.z = pid_control_V1_B.Gain;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_B.Vd1 = pid_control_V1_B.CL_h_OGE * pid_control_V1_B.hw;

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   * */
  pid_control_V1_B.BusAssignment.state.pose.orientation.w =
    pid_control_V1_B.alpha * pid_control_V1_B.beta * pid_control_V1_B.CL_w_OGE +
    pid_control_V1_B.Vd1 * pid_control_V1_B.CL_h_IGE;
  pid_control_V1_B.BusAssignment.state.pose.orientation.z = pid_control_V1_B.Vd1
    * pid_control_V1_B.CL_w_OGE - pid_control_V1_B.CL_h_IGE *
    pid_control_V1_B.alpha * pid_control_V1_B.beta;
  pid_control_V1_B.BusAssignment.state.pose.orientation.y =
    pid_control_V1_B.CL_h_OGE * pid_control_V1_B.CL_h_IGE *
    pid_control_V1_B.beta + pid_control_V1_B.hw * pid_control_V1_B.alpha *
    pid_control_V1_B.CL_w_OGE;
  pid_control_V1_B.BusAssignment.state.pose.orientation.x = pid_control_V1_B.hw *
    pid_control_V1_B.CL_h_IGE * pid_control_V1_B.alpha -
    pid_control_V1_B.CL_h_OGE * pid_control_V1_B.beta *
    pid_control_V1_B.CL_w_OGE;
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
  serverAvailableOnTime = ServCall_pid_control_V1_326.waitForServer(5.0);
  if (serverAvailableOnTime) {
    ServCall_pid_control_V1_326.call(&pid_control_V1_B.BusAssignment, &tmp);
  }

  /* End of MATLABSystem: '<S2>/ServiceCaller' */
  /* End of Outputs for SubSystem: '<Root>/Call Service' */

  /* Product: '<S12>/Product2' incorporates:
   *  Math: '<S12>/Square'
   *  Math: '<S12>/Square1'
   *  Math: '<S12>/Square2'
   *  Sqrt: '<S12>/Sqrt'
   *  Sum: '<S12>/Sum2'
   */
  pid_control_V1_B.Power = sqrt((pid_control_V1_B.x[0] * pid_control_V1_B.x[0] +
    pid_control_V1_B.x[1] * pid_control_V1_B.x[1]) + pid_control_V1_B.x[2] *
    pid_control_V1_B.x[2]) * pid_control_V1_B.XDOT[34];

  /* Gain: '<S12>/Gain3' */
  pid_control_V1_B.Gain3 = 0.001 * pid_control_V1_B.Power;
  if (tmp_0) {
  }

  /* Product: '<S12>/Divide' incorporates:
   *  Constant: '<S12>/thrust efficiency Cp?'
   */
  pid_control_V1_B.powerdemand = pid_control_V1_B.Gain3 / 0.248;
  if (tmp_0) {
  }

  /* Product: '<S12>/Divide1' */
  pid_control_V1_B.loadtorque = pid_control_V1_B.powerdemand /
    pid_control_V1_ConstB.motorspeed;
  if (tmp_0) {
  }

  /* Gain: '<S145>/ZeroGain' */
  pid_control_V1_B.alpha = 0.0 * pid_control_V1_B.SignPreSat_o;

  /* DeadZone: '<S147>/DeadZone' */
  if (pid_control_V1_B.SignPreSat_o > 1.0) {
    pid_control_V1_B.SignPreSat_o--;
  } else if (pid_control_V1_B.SignPreSat_o >= 0.0) {
    pid_control_V1_B.SignPreSat_o = 0.0;
  }

  /* End of DeadZone: '<S147>/DeadZone' */

  /* Gain: '<S152>/Integral Gain' */
  pid_control_V1_B.Switch_c *= 0.02;

  /* Signum: '<S145>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat_o)) {
    /* DataTypeConversion: '<S145>/DataTypeConv1' */
    i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat_o < 0.0) {
      /* DataTypeConversion: '<S145>/DataTypeConv1' */
      pid_control_V1_B.Vd1 = -1.0;
    } else {
      /* DataTypeConversion: '<S145>/DataTypeConv1' */
      pid_control_V1_B.Vd1 = (pid_control_V1_B.SignPreSat_o > 0.0);
    }

    /* DataTypeConversion: '<S145>/DataTypeConv1' */
    i = static_cast<int32_T>(fmod(pid_control_V1_B.Vd1, 256.0));
  }

  /* End of Signum: '<S145>/SignPreSat' */

  /* Signum: '<S145>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch_c)) {
    /* DataTypeConversion: '<S145>/DataTypeConv2' */
    tmp_1 = 0;
  } else {
    if (pid_control_V1_B.Switch_c < 0.0) {
      /* DataTypeConversion: '<S145>/DataTypeConv2' */
      pid_control_V1_B.Vd1 = -1.0;
    } else {
      /* DataTypeConversion: '<S145>/DataTypeConv2' */
      pid_control_V1_B.Vd1 = (pid_control_V1_B.Switch_c > 0.0);
    }

    /* DataTypeConversion: '<S145>/DataTypeConv2' */
    tmp_1 = static_cast<int32_T>(fmod(pid_control_V1_B.Vd1, 256.0));
  }

  /* End of Signum: '<S145>/SignPreIntegrator' */

  /* DataTypeConversion: '<S145>/DataTypeConv1' */
  if (i < 0) {
    i = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(i))));
  }

  /* DataTypeConversion: '<S145>/DataTypeConv2' */
  if (tmp_1 < 0) {
    tmp_1 = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(tmp_1))));
  }

  /* Logic: '<S145>/AND3' incorporates:
   *  DataTypeConversion: '<S145>/DataTypeConv1'
   *  DataTypeConversion: '<S145>/DataTypeConv2'
   *  RelationalOperator: '<S145>/Equal1'
   *  RelationalOperator: '<S145>/NotEqual'
   */
  pid_control_V1_B.AND3 = ((pid_control_V1_B.alpha !=
    pid_control_V1_B.SignPreSat_o) && (i == tmp_1));
  if (tmp_0) {
    /* Memory: '<S145>/Memory' */
    pid_control_V1_B.Memory = pid_control_V1_DW.Memory_PreviousInput;
  }

  /* Switch: '<S145>/Switch' */
  if (pid_control_V1_B.Memory) {
    /* Gain: '<S152>/Integral Gain' incorporates:
     *  Constant: '<S145>/Constant1'
     *  Switch: '<S145>/Switch'
     */
    pid_control_V1_B.Switch_c = 0.0;
  }

  /* End of Switch: '<S145>/Switch' */

  /* Gain: '<S44>/Integral Gain' incorporates:
   *  Constant: '<Root>/Constant4'
   *  Sum: '<Root>/Sum4'
   */
  pid_control_V1_B.IntegralGain = -(0.0 - pid_control_V1_B.x[6]);

  /* Gain: '<S204>/Integral Gain' */
  pid_control_V1_B.IntegralGain_a = -0.5 * pid_control_V1_B.Sum5;

  /* Gain: '<S91>/ZeroGain' */
  pid_control_V1_B.Sum5 = 0.0 * pid_control_V1_B.SignPreSat;

  /* DeadZone: '<S93>/DeadZone' */
  if (pid_control_V1_B.SignPreSat > 0.3490658503988659) {
    pid_control_V1_B.SignPreSat -= 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat >= -0.3490658503988659) {
    pid_control_V1_B.SignPreSat = 0.0;
  } else {
    pid_control_V1_B.SignPreSat -= -0.3490658503988659;
  }

  /* End of DeadZone: '<S93>/DeadZone' */

  /* Gain: '<S98>/Integral Gain' */
  pid_control_V1_B.Switch_a *= -5.0;

  /* Signum: '<S91>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat)) {
    /* DataTypeConversion: '<S91>/DataTypeConv1' */
    i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat < 0.0) {
      /* DataTypeConversion: '<S91>/DataTypeConv1' */
      pid_control_V1_B.Vd1 = -1.0;
    } else {
      /* DataTypeConversion: '<S91>/DataTypeConv1' */
      pid_control_V1_B.Vd1 = (pid_control_V1_B.SignPreSat > 0.0);
    }

    /* DataTypeConversion: '<S91>/DataTypeConv1' */
    i = static_cast<int32_T>(fmod(pid_control_V1_B.Vd1, 256.0));
  }

  /* End of Signum: '<S91>/SignPreSat' */

  /* Signum: '<S91>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch_a)) {
    /* DataTypeConversion: '<S91>/DataTypeConv2' */
    tmp_1 = 0;
  } else {
    if (pid_control_V1_B.Switch_a < 0.0) {
      /* DataTypeConversion: '<S91>/DataTypeConv2' */
      pid_control_V1_B.Vd1 = -1.0;
    } else {
      /* DataTypeConversion: '<S91>/DataTypeConv2' */
      pid_control_V1_B.Vd1 = (pid_control_V1_B.Switch_a > 0.0);
    }

    /* DataTypeConversion: '<S91>/DataTypeConv2' */
    tmp_1 = static_cast<int32_T>(fmod(pid_control_V1_B.Vd1, 256.0));
  }

  /* End of Signum: '<S91>/SignPreIntegrator' */

  /* DataTypeConversion: '<S91>/DataTypeConv1' */
  if (i < 0) {
    i = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(i))));
  }

  /* DataTypeConversion: '<S91>/DataTypeConv2' */
  if (tmp_1 < 0) {
    tmp_1 = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(tmp_1))));
  }

  /* Logic: '<S91>/AND3' incorporates:
   *  DataTypeConversion: '<S91>/DataTypeConv1'
   *  DataTypeConversion: '<S91>/DataTypeConv2'
   *  RelationalOperator: '<S91>/Equal1'
   *  RelationalOperator: '<S91>/NotEqual'
   */
  pid_control_V1_B.AND3_j = ((pid_control_V1_B.Sum5 !=
    pid_control_V1_B.SignPreSat) && (i == tmp_1));
  if (tmp_0) {
    /* Memory: '<S91>/Memory' */
    pid_control_V1_B.Memory_n = pid_control_V1_DW.Memory_PreviousInput_d;
  }

  /* Switch: '<S91>/Switch' */
  if (pid_control_V1_B.Memory_n) {
    /* Gain: '<S98>/Integral Gain' incorporates:
     *  Constant: '<S91>/Constant1'
     *  Switch: '<S91>/Switch'
     */
    pid_control_V1_B.Switch_a = 0.0;
  }

  /* End of Switch: '<S91>/Switch' */

  /* Gain: '<S12>/Gain1' incorporates:
   *  Integrator: '<S12>/Integrator1'
   */
  pid_control_V1_B.EnergykWh = 2.7777777777777776E-7 *
    pid_control_V1_X.Integrator1_CSTATE;
  if (tmp_0) {
  }

  if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
    if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
      /* Update for UnitDelay: '<Root>/Unit Delay1' */
      pid_control_V1_DW.UnitDelay1_DSTATE = pid_control_V1_B.Switch1;

      /* Update for UnitDelay: '<Root>/Unit Delay2' */
      pid_control_V1_DW.UnitDelay2_DSTATE = pid_control_V1_B.Switch2;

      /* Update for UnitDelay: '<Root>/Unit Delay' */
      pid_control_V1_DW.UnitDelay_DSTATE = pid_control_V1_B.Switch;

      /* Update for Memory: '<S145>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput = pid_control_V1_B.AND3;

      /* Update for Memory: '<S91>/Memory' */
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

  /* Derivatives for Integrator: '<S12>/Integrator' */
  memcpy(&_rtXdot->Integrator_CSTATE[0], &pid_control_V1_B.XDOT[0], 12U * sizeof
         (real_T));

  /* Derivatives for Integrator: '<S101>/Integrator' */
  _rtXdot->Integrator_CSTATE_p = pid_control_V1_B.Switch_a;

  /* Derivatives for Integrator: '<S96>/Filter' */
  _rtXdot->Filter_CSTATE = pid_control_V1_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S207>/Integrator' */
  _rtXdot->Integrator_CSTATE_d = pid_control_V1_B.IntegralGain_a;

  /* Derivatives for Integrator: '<S202>/Filter' */
  _rtXdot->Filter_CSTATE_f = pid_control_V1_B.FilterCoefficient_p;

  /* Derivatives for Integrator: '<S47>/Integrator' */
  _rtXdot->Integrator_CSTATE_m = pid_control_V1_B.IntegralGain;

  /* Derivatives for Integrator: '<S42>/Filter' */
  _rtXdot->Filter_CSTATE_g = pid_control_V1_B.FilterCoefficient_c;

  /* Derivatives for Integrator: '<S155>/Integrator' */
  _rtXdot->Integrator_CSTATE_b = pid_control_V1_B.Switch_c;

  /* Derivatives for Integrator: '<S150>/Filter' */
  _rtXdot->Filter_CSTATE_fy = pid_control_V1_B.FilterCoefficient_g;

  /* Derivatives for Integrator: '<S12>/Integrator1' */
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

  /* Start for MATLABSystem: '<S10>/SourceBlock' */
  pid_control_V1_DW.obj_j.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_j.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_a = true;
  pid_control_V1_DW.obj_j.isSetupComplete = false;
  pid_control_V1_DW.obj_j.isInitialized = 1;
  pid_cont_Subscriber_setupImpl_o(&pid_control_V1_DW.obj_j);
  pid_control_V1_DW.obj_j.isSetupComplete = true;

  /* Start for MATLABSystem: '<S11>/SourceBlock' */
  pid_control_V1_DW.obj_k.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_k.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty = true;
  pid_control_V1_DW.obj_k.isSetupComplete = false;
  pid_control_V1_DW.obj_k.isInitialized = 1;
  pid_con_Subscriber_setupImpl_on(&pid_control_V1_DW.obj_k);
  pid_control_V1_DW.obj_k.isSetupComplete = true;

  /* Start for MATLABSystem: '<S9>/SourceBlock' */
  pid_control_V1_DW.obj_h.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_h.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_n = true;
  pid_control_V1_DW.obj_h.isSetupComplete = false;
  pid_control_V1_DW.obj_h.isInitialized = 1;
  pid_contro_Subscriber_setupImpl(&pid_control_V1_DW.obj_h);
  pid_control_V1_DW.obj_h.isSetupComplete = true;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_DW.objisempty_d = true;
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

  /* InitializeConditions for Integrator: '<S12>/Integrator' */
  memcpy(&pid_control_V1_X.Integrator_CSTATE[0],
         &pid_control_V1_ConstP.Integrator_IC[0], 12U * sizeof(real_T));

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay1' */
  pid_control_V1_DW.UnitDelay1_DSTATE = 0.0087266462599716477;

  /* InitializeConditions for Integrator: '<S101>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S96>/Filter' */
  pid_control_V1_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay2' */
  pid_control_V1_DW.UnitDelay2_DSTATE = 0.26179938779914941;

  /* InitializeConditions for Integrator: '<S207>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_d = 0.0;

  /* InitializeConditions for Integrator: '<S202>/Filter' */
  pid_control_V1_X.Filter_CSTATE_f = 0.0;

  /* InitializeConditions for Integrator: '<S47>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_m = 0.0;

  /* InitializeConditions for Integrator: '<S42>/Filter' */
  pid_control_V1_X.Filter_CSTATE_g = 0.0;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay' */
  pid_control_V1_DW.UnitDelay_DSTATE = 1.0;

  /* InitializeConditions for Integrator: '<S155>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_b = 0.3;

  /* InitializeConditions for Integrator: '<S150>/Filter' */
  pid_control_V1_X.Filter_CSTATE_fy = 0.0;

  /* InitializeConditions for Integrator: '<S12>/Integrator1' */
  pid_control_V1_X.Integrator1_CSTATE = 0.0;

  /* SystemInitialize for Enabled SubSystem: '<S10>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_c);

  /* End of SystemInitialize for SubSystem: '<S10>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S11>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_a);

  /* End of SystemInitialize for SubSystem: '<S11>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S9>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem);

  /* End of SystemInitialize for SubSystem: '<S9>/Enabled Subsystem' */
}

/* Model terminate function */
void pid_control_V1::terminate()
{
  /* Terminate for MATLABSystem: '<S10>/SourceBlock' */
  if (!pid_control_V1_DW.obj_j.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_j.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_j.isInitialized == 1) &&
        pid_control_V1_DW.obj_j.isSetupComplete) {
      Sub_pid_control_V1_370.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S10>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S11>/SourceBlock' */
  if (!pid_control_V1_DW.obj_k.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_k.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_k.isInitialized == 1) &&
        pid_control_V1_DW.obj_k.isSetupComplete) {
      Sub_pid_control_V1_377.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S11>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S9>/SourceBlock' */
  if (!pid_control_V1_DW.obj_h.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_h.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_h.isInitialized == 1) &&
        pid_control_V1_DW.obj_h.isSetupComplete) {
      Sub_pid_control_V1_366.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S9>/SourceBlock' */
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
