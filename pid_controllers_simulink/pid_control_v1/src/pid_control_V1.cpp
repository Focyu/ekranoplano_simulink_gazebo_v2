/*
 * pid_control_V1.cpp
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.96
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Mon Mar 23 00:33:26 2026
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
#include "pid_control_V1_private.h"
#include <emmintrin.h>
#include <math.h>

extern "C"
{

#include "rt_nonfinite.h"

}

#include "rmw/qos_profiles.h"
#include <stddef.h>
#include "rt_defines.h"

uint32_T plook_bincpa(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                      *fraction, uint32_T *prevIndex)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'on'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = 0.0;
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d_prevIdx(u, bp, *prevIndex, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex;
    *fraction = 0.0;
  }

  *prevIndex = bpIndex;
  return bpIndex;
}

real_T intrp2d_la_pw(const uint32_T bpIndex[], const real_T frac[], const real_T
                     table[], const uint32_T stride, const uint32_T maxIndex[])
{
  real_T y;
  real_T yL_0d0;
  uint32_T offset_1d;

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'on'
     Overflow mode: 'portable wrapping'
   */
  offset_1d = bpIndex[1U] * stride + bpIndex[0U];
  if (bpIndex[0U] == maxIndex[0U]) {
    y = table[offset_1d];
  } else {
    yL_0d0 = table[offset_1d];
    y = (table[offset_1d + 1U] - yL_0d0) * frac[0U] + yL_0d0;
  }

  if (bpIndex[1U] == maxIndex[1U]) {
  } else {
    offset_1d += stride;
    if (bpIndex[0U] == maxIndex[0U]) {
      yL_0d0 = table[offset_1d];
    } else {
      yL_0d0 = table[offset_1d];
      yL_0d0 += (table[offset_1d + 1U] - yL_0d0) * frac[0U];
    }

    y += (yL_0d0 - y) * frac[1U];
  }

  return y;
}

uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIndex;
  uint32_T found;
  uint32_T iLeft;
  uint32_T iRght;

  /* Binary Search using Previous Index */
  bpIndex = startIndex;
  iLeft = 0U;
  iRght = maxIndex;
  found = 0U;
  while (found == 0U) {
    if (u < bp[bpIndex]) {
      iRght = bpIndex - 1U;
      bpIndex = ((bpIndex + iLeft) - 1U) >> 1U;
    } else if (u < bp[bpIndex + 1U]) {
      found = 1U;
    } else {
      iLeft = bpIndex + 1U;
      bpIndex = ((bpIndex + iRght) + 1U) >> 1U;
    }
  }

  return bpIndex;
}

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
  int_T nXc = 42;
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
 *    '<S11>/Enabled Subsystem'
 *    '<S12>/Enabled Subsystem'
 *    '<S13>/Enabled Subsystem'
 *    '<S293>/Enabled Subsystem'
 *    '<S294>/Enabled Subsystem'
 */
void pid_control_V1::pid_contr_EnabledSubsystem_Init
  (B_EnabledSubsystem_pid_contro_T *localB)
{
  /* SystemInitialize for SignalConversion generated from: '<S283>/In1' */
  memset(&localB->In1, 0, sizeof(SL_Bus_std_msgs_Float64));
}

/*
 * Output and update for enable system:
 *    '<S11>/Enabled Subsystem'
 *    '<S12>/Enabled Subsystem'
 *    '<S13>/Enabled Subsystem'
 *    '<S293>/Enabled Subsystem'
 *    '<S294>/Enabled Subsystem'
 */
void pid_control_V1::pid_control_V1_EnabledSubsystem(boolean_T rtu_Enable, const
  SL_Bus_std_msgs_Float64 *rtu_In1, B_EnabledSubsystem_pid_contro_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S11>/Enabled Subsystem' incorporates:
   *  EnablePort: '<S283>/Enable'
   */
  if (rtu_Enable) {
    /* SignalConversion generated from: '<S283>/In1' */
    localB->In1 = *rtu_In1;
  }

  /* End of Outputs for SubSystem: '<S11>/Enabled Subsystem' */
}

/*
 * System initialize for enable system:
 *    '<S14>/Enabled Subsystem'
 *    '<S15>/Enabled Subsystem'
 */
void pid_control_V1::pid_con_EnabledSubsystem_d_Init
  (B_EnabledSubsystem_pid_cont_p_T *localB)
{
  /* SystemInitialize for SignalConversion generated from: '<S286>/In1' */
  memset(&localB->In1, 0, sizeof(SL_Bus_std_msgs_Int64));
}

/*
 * Output and update for enable system:
 *    '<S14>/Enabled Subsystem'
 *    '<S15>/Enabled Subsystem'
 */
void pid_control_V1::pid_control__EnabledSubsystem_h(boolean_T rtu_Enable, const
  SL_Bus_std_msgs_Int64 *rtu_In1, B_EnabledSubsystem_pid_cont_p_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S14>/Enabled Subsystem' incorporates:
   *  EnablePort: '<S286>/Enable'
   */
  if (rtu_Enable) {
    /* SignalConversion generated from: '<S286>/In1' */
    localB->In1 = *rtu_In1;
  }

  /* End of Outputs for SubSystem: '<S14>/Enabled Subsystem' */
}

/*
 * System initialize for enable system:
 *    '<S291>/Enabled Subsystem'
 *    '<S292>/Enabled Subsystem'
 */
void pid_control_V1::pid_con_EnabledSubsystem_i_Init
  (B_EnabledSubsystem_pid_cont_n_T *localB)
{
  /* SystemInitialize for SignalConversion generated from: '<S334>/In1' */
  memset(&localB->In1, 0, sizeof(SL_Bus_std_msgs_Bool));
}

/*
 * Output and update for enable system:
 *    '<S291>/Enabled Subsystem'
 *    '<S292>/Enabled Subsystem'
 */
void pid_control_V1::pid_control__EnabledSubsystem_p(boolean_T rtu_Enable, const
  SL_Bus_std_msgs_Bool *rtu_In1, B_EnabledSubsystem_pid_cont_n_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S291>/Enabled Subsystem' incorporates:
   *  EnablePort: '<S334>/Enable'
   */
  if (rtu_Enable) {
    /* SignalConversion generated from: '<S334>/In1' */
    localB->In1 = *rtu_In1;
  }

  /* End of Outputs for SubSystem: '<S291>/Enabled Subsystem' */
}

void pid_control_V1::pid_con_Subscriber_setupImpl_on(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[17] = "/setpoint/altura";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S13>/SourceBlock' */
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
    /* Start for MATLABSystem: '<S13>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_cv[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_435.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_cv
    [0], qos_profile);
}

void pid_control_V1::pid_contro_Subscriber_setupImpl(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[21] = "/setpoint/waypoint_x";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S11>/SourceBlock' */
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
  for (int32_T i = 0; i < 21; i++) {
    /* Start for MATLABSystem: '<S11>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_b[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_466.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_b[0],
    qos_profile);
}

void pid_control_V1::pid_co_Subscriber_setupImpl_onh(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[21] = "/setpoint/waypoint_y";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S14>/SourceBlock' */
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
  for (int32_T i = 0; i < 21; i++) {
    /* Start for MATLABSystem: '<S14>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_cx[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_467.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_cx
    [0], qos_profile);
}

void pid_control_V1::pid_c_Subscriber_setupImpl_onhg(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[15];
  static const char_T b_zeroDelimTopic_0[15] = "/setpoint/mode";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S15>/SourceBlock' */
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
  for (int32_T i = 0; i < 15; i++) {
    /* Start for MATLABSystem: '<S15>/SourceBlock' */
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_pid_control_V1_476.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void pid_control_V1::pid_cont_Subscriber_setupImpl_o(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[14];
  static const char_T b_zeroDelimTopic_0[14] = "/setpoint/yaw";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S12>/SourceBlock' */
  pid_control_V1_B.deadline_f.sec = 0.0;
  pid_control_V1_B.deadline_f.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)10.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE,
                 pid_control_V1_B.deadline_f, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 14; i++) {
    /* Start for MATLABSystem: '<S12>/SourceBlock' */
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_pid_control_V1_377.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
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

void pid_control_V1::pid__Subscriber_setupImpl_onhgd(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[22] = "/setpoint/turbulencia";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S291>/SourceBlock' */
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
  for (int32_T i = 0; i < 22; i++) {
    /* Start for MATLABSystem: '<S291>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_k[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_417.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_k[0],
    qos_profile);
}

void pid_control_V1::pid_Subscriber_setupImpl_onhgd0(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[22] = "/setpoint/turbulencia";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S292>/SourceBlock' */
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
  for (int32_T i = 0; i < 22; i++) {
    /* Start for MATLABSystem: '<S292>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_c[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_423.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_c[0],
    qos_profile);
}

void pid_control_V1::pi_Subscriber_setupImpl_onhgd03(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[12];
  static const char_T b_zeroDelimTopic_0[12] = "/olas/heave";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S293>/SourceBlock' */
  pid_control_V1_B.deadline_g.sec = 0.0;
  pid_control_V1_B.deadline_g.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)10.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_RELIABLE,
                 pid_control_V1_B.deadline_g, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 12; i++) {
    /* Start for MATLABSystem: '<S293>/SourceBlock' */
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_pid_control_V1_443.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void pid_control_V1::p_Subscriber_setupImpl_onhgd03r(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  static const char_T b_zeroDelimTopic[17] = "/olas/pitch_rate";
  qos_profile = rmw_qos_profile_default;

  /* Start for MATLABSystem: '<S294>/SourceBlock' */
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
    /* Start for MATLABSystem: '<S294>/SourceBlock' */
    pid_control_V1_B.b_zeroDelimTopic_p[i] = b_zeroDelimTopic[i];
  }

  Sub_pid_control_V1_445.createSubscriber(&pid_control_V1_B.b_zeroDelimTopic_p[0],
    qos_profile);
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T hi;
  uint32_T lo;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return static_cast<real_T>(*u) * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T si;
  real_T sr;
  real_T y;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
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
  /* local block i/o variables */
  SL_Bus_std_msgs_Float64 rtb_SourceBlock_o2_d;
  SL_Bus_std_msgs_Float64 rtb_SourceBlock_o2_l;
  SL_Bus_std_msgs_Bool rtb_SourceBlock_o2_dd;
  SL_Bus_std_msgs_Bool rtb_SourceBlock_o2_j;
  __m128d tmp_2;
  SL_Bus_gazebo_msgs_SetEntityStateResponse tmp;
  int32_T i;
  int32_T tmp_1;
  int8_T rtAction;
  int8_T rtPrevAction;
  boolean_T serverAvailableOnTime;
  boolean_T tmp_0;
  static const uint8_T b[11] = { 101U, 107U, 114U, 97U, 110U, 111U, 112U, 108U,
    97U, 110U, 111U };

  static const uint8_T b_0[5] = { 119U, 111U, 114U, 108U, 100U };

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

  /* Outputs for Enabled SubSystem: '<S297>/Hugw(s)' incorporates:
   *  EnablePort: '<S310>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S296>/Hrgw' incorporates:
   *  EnablePort: '<S309>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S297>/Hvgw(s)' incorporates:
   *  EnablePort: '<S311>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S296>/Hqgw' incorporates:
   *  EnablePort: '<S308>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S297>/Hwgw(s)' incorporates:
   *  EnablePort: '<S312>/Enable'
   */
  /* Outputs for Enabled SubSystem: '<S296>/Hpgw' incorporates:
   *  EnablePort: '<S307>/Enable'
   */
  tmp_0 = rtmIsMajorTimeStep((&pid_control_V1_M));

  /* End of Outputs for SubSystem: '<S296>/Hpgw' */
  /* End of Outputs for SubSystem: '<S297>/Hwgw(s)' */
  /* End of Outputs for SubSystem: '<S296>/Hqgw' */
  /* End of Outputs for SubSystem: '<S297>/Hvgw(s)' */
  /* End of Outputs for SubSystem: '<S296>/Hrgw' */
  /* End of Outputs for SubSystem: '<S297>/Hugw(s)' */
  if (tmp_0) {
    /* MATLABSystem: '<S13>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_o = Sub_pid_control_V1_435.getLatestMessage(
      &rtb_SourceBlock_o2_d);

    /* Outputs for Enabled SubSystem: '<S13>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_o,
      &rtb_SourceBlock_o2_d, &pid_control_V1_B.EnabledSubsystem_b);

    /* End of Outputs for SubSystem: '<S13>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch3' */
    if (pid_control_V1_B.SourceBlock_o1_o) {
      /* Switch: '<Root>/Switch3' */
      pid_control_V1_B.Switch3 = pid_control_V1_B.EnabledSubsystem_b.In1.data;
    } else {
      /* Switch: '<Root>/Switch3' incorporates:
       *  UnitDelay: '<Root>/Unit Delay3'
       */
      pid_control_V1_B.Switch3 = pid_control_V1_DW.UnitDelay3_DSTATE;
    }

    /* End of Switch: '<Root>/Switch3' */
  }

  /* Integrator: '<S16>/Integrator' */
  memcpy(&pid_control_V1_B.x[0], &pid_control_V1_X.Integrator_CSTATE[0], 12U *
         sizeof(real_T));

  /* Gain: '<Root>/Gain' */
  pid_control_V1_B.Gain = -pid_control_V1_B.x[11];

  /* Gain: '<S104>/Integral Gain' incorporates:
   *  Gain: '<Root>/Gain4'
   *  Sum: '<Root>/Sum2'
   */
  pid_control_V1_B.Switch_k = pid_control_V1_B.Switch3 - (-pid_control_V1_B.x[11]);

  /* Gain: '<S110>/Filter Coefficient' incorporates:
   *  Gain: '<S100>/Derivative Gain'
   *  Integrator: '<S102>/Filter'
   *  Sum: '<S102>/SumD'
   */
  pid_control_V1_B.FilterCoefficient = (0.5 * pid_control_V1_B.Switch_k -
    pid_control_V1_X.Filter_CSTATE) * 100.0;

  /* Sum: '<S116>/Sum' incorporates:
   *  Integrator: '<S107>/Integrator'
   */
  pid_control_V1_B.SignPreSat = (pid_control_V1_B.Switch_k +
    pid_control_V1_X.Integrator_CSTATE_n) + pid_control_V1_B.FilterCoefficient;

  /* Saturate: '<S114>/Saturation' */
  if (pid_control_V1_B.SignPreSat > 20.0) {
    /* Saturate: '<S114>/Saturation' */
    pid_control_V1_B.Saturation = 20.0;
  } else if (pid_control_V1_B.SignPreSat < 0.0) {
    /* Saturate: '<S114>/Saturation' */
    pid_control_V1_B.Saturation = 0.0;
  } else {
    /* Saturate: '<S114>/Saturation' */
    pid_control_V1_B.Saturation = pid_control_V1_B.SignPreSat;
  }

  /* End of Saturate: '<S114>/Saturation' */

  /* Gain: '<S56>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S46>/Derivative Gain'
   *  Integrator: '<S48>/Filter'
   *  Sum: '<Root>/Sum4'
   *  Sum: '<S48>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_c = ((0.0 - pid_control_V1_B.x[6]) * -0.1 -
    pid_control_V1_X.Filter_CSTATE_g) * 100.0;

  /* Sum: '<S62>/Sum' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Gain: '<S58>/Proportional Gain'
   *  Integrator: '<S53>/Integrator'
   *  Sum: '<Root>/Sum4'
   */
  pid_control_V1_B.SignPreSat_a = (-(0.0 - pid_control_V1_B.x[6]) +
    pid_control_V1_X.Integrator_CSTATE_m) + pid_control_V1_B.FilterCoefficient_c;

  /* Saturate: '<S60>/Saturation' */
  if (pid_control_V1_B.SignPreSat_a > 0.3490658503988659) {
    /* Saturate: '<S60>/Saturation' */
    pid_control_V1_B.Saturation_k = 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat_a < -0.3490658503988659) {
    /* Saturate: '<S60>/Saturation' */
    pid_control_V1_B.Saturation_k = -0.3490658503988659;
  } else {
    /* Saturate: '<S60>/Saturation' */
    pid_control_V1_B.Saturation_k = pid_control_V1_B.SignPreSat_a;
  }

  /* End of Saturate: '<S60>/Saturation' */

  /* Saturate: '<Root>/Saturation' */
  if (pid_control_V1_B.Saturation > 0.13962634015954636) {
    /* Saturate: '<Root>/Saturation' */
    pid_control_V1_B.Saturation_i = 0.13962634015954636;
  } else if (pid_control_V1_B.Saturation < -0.034906585039886591) {
    /* Saturate: '<Root>/Saturation' */
    pid_control_V1_B.Saturation_i = -0.034906585039886591;
  } else {
    /* Saturate: '<Root>/Saturation' */
    pid_control_V1_B.Saturation_i = pid_control_V1_B.Saturation;
  }

  /* End of Saturate: '<Root>/Saturation' */

  /* Sum: '<Root>/Sum1' */
  pid_control_V1_B.Sum1_g = pid_control_V1_B.Saturation_i - pid_control_V1_B.x[7];

  /* Gain: '<S162>/Filter Coefficient' incorporates:
   *  Gain: '<S152>/Derivative Gain'
   *  Integrator: '<S154>/Filter'
   *  Sum: '<S154>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_m = (-pid_control_V1_B.Sum1_g -
    pid_control_V1_X.Filter_CSTATE_m) * 100.0;

  /* Sum: '<S168>/Sum' incorporates:
   *  Gain: '<S164>/Proportional Gain'
   *  Integrator: '<S159>/Integrator'
   */
  pid_control_V1_B.Sum_hl = (-pid_control_V1_B.Sum1_g +
    pid_control_V1_X.Integrator_CSTATE_p) + pid_control_V1_B.FilterCoefficient_m;

  /* Saturate: '<S166>/Saturation' */
  if (pid_control_V1_B.Sum_hl > 0.3490658503988659) {
    /* Saturate: '<S166>/Saturation' */
    pid_control_V1_B.Saturation_f = 0.3490658503988659;
  } else if (pid_control_V1_B.Sum_hl < -0.3490658503988659) {
    /* Saturate: '<S166>/Saturation' */
    pid_control_V1_B.Saturation_f = -0.3490658503988659;
  } else {
    /* Saturate: '<S166>/Saturation' */
    pid_control_V1_B.Saturation_f = pid_control_V1_B.Sum_hl;
  }

  /* End of Saturate: '<S166>/Saturation' */
  if (tmp_0) {
    /* MATLABSystem: '<S11>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_l = Sub_pid_control_V1_466.getLatestMessage(
      &rtb_SourceBlock_o2_l);

    /* Outputs for Enabled SubSystem: '<S11>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_l,
      &rtb_SourceBlock_o2_l, &pid_control_V1_B.EnabledSubsystem);

    /* End of Outputs for SubSystem: '<S11>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch4' */
    if (pid_control_V1_B.SourceBlock_o1_l) {
      /* Switch: '<Root>/Switch4' */
      pid_control_V1_B.Switch4 = pid_control_V1_B.EnabledSubsystem.In1.data;
    } else {
      /* Switch: '<Root>/Switch4' incorporates:
       *  UnitDelay: '<Root>/Unit Delay4'
       */
      pid_control_V1_B.Switch4 = pid_control_V1_DW.UnitDelay4_DSTATE;
    }

    /* End of Switch: '<Root>/Switch4' */

    /* MATLABSystem: '<S14>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_g = Sub_pid_control_V1_467.getLatestMessage(
      &pid_control_V1_B.SourceBlock_o2_g);

    /* Outputs for Enabled SubSystem: '<S14>/Enabled Subsystem' */
    pid_control__EnabledSubsystem_h(pid_control_V1_B.SourceBlock_o1_g,
      &pid_control_V1_B.SourceBlock_o2_g, &pid_control_V1_B.EnabledSubsystem_h);

    /* End of Outputs for SubSystem: '<S14>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch5' */
    if (pid_control_V1_B.SourceBlock_o1_g) {
      /* Switch: '<Root>/Switch5' */
      pid_control_V1_B.Switch5 = pid_control_V1_B.EnabledSubsystem_h.In1.data;
    } else {
      /* Switch: '<Root>/Switch5' incorporates:
       *  UnitDelay: '<Root>/Unit Delay5'
       */
      pid_control_V1_B.Switch5 = pid_control_V1_DW.UnitDelay5_DSTATE;
    }

    /* End of Switch: '<Root>/Switch5' */
  }

  /* MATLAB Function: '<Root>/MATLAB Function2' */
  tmp_2 = _mm_sub_pd(_mm_set_pd(pid_control_V1_B.Switch5,
    pid_control_V1_B.Switch4), _mm_loadu_pd(&pid_control_V1_B.x[9]));
  _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);

  /* MATLAB Function: '<Root>/MATLAB Function2' */
  if ((fabs(pid_control_V1_B.dv[0]) < 0.01) && (fabs(pid_control_V1_B.dv[1]) <
       0.01)) {
    pid_control_V1_B.Switch = pid_control_V1_B.x[8];
  } else {
    pid_control_V1_B.chi = rt_atan2d_snf(pid_control_V1_B.dv[1],
      pid_control_V1_B.dv[0]);
    pid_control_V1_B.Switch = rt_atan2d_snf(sin(pid_control_V1_B.chi -
      pid_control_V1_B.x[8]), cos(pid_control_V1_B.chi - pid_control_V1_B.x[8]))
      + pid_control_V1_B.x[8];
  }

  if (tmp_0) {
    /* MATLABSystem: '<S15>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_d = Sub_pid_control_V1_476.getLatestMessage(
      &pid_control_V1_B.SourceBlock_o2);

    /* Outputs for Enabled SubSystem: '<S15>/Enabled Subsystem' */
    pid_control__EnabledSubsystem_h(pid_control_V1_B.SourceBlock_o1_d,
      &pid_control_V1_B.SourceBlock_o2, &pid_control_V1_B.EnabledSubsystem_bk);

    /* End of Outputs for SubSystem: '<S15>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch6' */
    if (pid_control_V1_B.SourceBlock_o1_d) {
      /* Switch: '<Root>/Switch6' */
      pid_control_V1_B.Switch6 = pid_control_V1_B.EnabledSubsystem_bk.In1.data;
    } else {
      /* Switch: '<Root>/Switch6' incorporates:
       *  UnitDelay: '<Root>/Unit Delay6'
       */
      pid_control_V1_B.Switch6 = pid_control_V1_DW.UnitDelay6_DSTATE;
    }

    /* End of Switch: '<Root>/Switch6' */

    /* MATLABSystem: '<S12>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_gx = Sub_pid_control_V1_377.getLatestMessage
      (&pid_control_V1_B.SourceBlock_o2_g2);

    /* Outputs for Enabled SubSystem: '<S12>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_gx,
      &pid_control_V1_B.SourceBlock_o2_g2, &pid_control_V1_B.EnabledSubsystem_a);

    /* End of Outputs for SubSystem: '<S12>/Enabled Subsystem' */

    /* Switch: '<Root>/Switch2' */
    if (pid_control_V1_B.SourceBlock_o1_gx) {
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

  /* Switch: '<Root>/Switch' */
  if (!(pid_control_V1_B.Switch6 > 0.0)) {
    /* Switch: '<Root>/Switch' */
    pid_control_V1_B.Switch = pid_control_V1_B.Switch2;
  }

  /* End of Switch: '<Root>/Switch' */

  /* Sum: '<Root>/Sum5' */
  pid_control_V1_B.Sum5 = pid_control_V1_B.Switch - pid_control_V1_B.x[8];

  /* Gain: '<S214>/Filter Coefficient' incorporates:
   *  Gain: '<S204>/Derivative Gain'
   *  Integrator: '<S206>/Filter'
   *  Sum: '<S206>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_p = (-pid_control_V1_B.Sum5 -
    pid_control_V1_X.Filter_CSTATE_f) * 100.0;

  /* Sum: '<S220>/Sum' incorporates:
   *  Gain: '<S216>/Proportional Gain'
   *  Integrator: '<S211>/Integrator'
   */
  pid_control_V1_B.Saturation_m = (-3.0 * pid_control_V1_B.Sum5 +
    pid_control_V1_X.Integrator_CSTATE_d) + pid_control_V1_B.FilterCoefficient_p;

  /* Saturate: '<S218>/Saturation' */
  if (pid_control_V1_B.Saturation_m > 0.26179938779914941) {
    /* Sum: '<S220>/Sum' incorporates:
     *  Saturate: '<S218>/Saturation'
     */
    pid_control_V1_B.Saturation_m = 0.26179938779914941;
  } else if (pid_control_V1_B.Saturation_m < -0.26179938779914941) {
    /* Sum: '<S220>/Sum' incorporates:
     *  Saturate: '<S218>/Saturation'
     */
    pid_control_V1_B.Saturation_m = -0.26179938779914941;
  }

  /* End of Saturate: '<S218>/Saturation' */
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
  pid_control_V1_B.sina = sin(pid_control_V1_B.wbe_b[0]);
  pid_control_V1_B.sinb = sin(pid_control_V1_B.wbe_b[1]);
  pid_control_V1_B.sinc = sin(pid_control_V1_B.wbe_b[2]);
  pid_control_V1_B.cosa = cos(pid_control_V1_B.wbe_b[0]);
  pid_control_V1_B.cosb = cos(pid_control_V1_B.wbe_b[1]);
  pid_control_V1_B.cosc = cos(pid_control_V1_B.wbe_b[2]);

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  Gain: '<Root>/Gain3'
   */
  pid_control_V1_B.BusAssignment.state.pose.position.x = pid_control_V1_B.x[9];
  pid_control_V1_B.BusAssignment.state.pose.position.y = -pid_control_V1_B.x[10];
  pid_control_V1_B.BusAssignment.state.pose.position.z = pid_control_V1_B.Gain;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_B.chi = pid_control_V1_B.cosa * pid_control_V1_B.cosb;

  /* BusAssignment: '<Root>/Bus Assignment' incorporates:
   *  MATLABSystem: '<Root>/Coordinate Transformation Conversion'
   * */
  pid_control_V1_B.BusAssignment.state.pose.orientation.w =
    pid_control_V1_B.sina * pid_control_V1_B.sinb * pid_control_V1_B.sinc +
    pid_control_V1_B.chi * pid_control_V1_B.cosc;
  pid_control_V1_B.BusAssignment.state.pose.orientation.z = pid_control_V1_B.chi
    * pid_control_V1_B.sinc - pid_control_V1_B.cosc * pid_control_V1_B.sina *
    pid_control_V1_B.sinb;
  pid_control_V1_B.BusAssignment.state.pose.orientation.y =
    pid_control_V1_B.cosa * pid_control_V1_B.cosc * pid_control_V1_B.sinb +
    pid_control_V1_B.cosb * pid_control_V1_B.sina * pid_control_V1_B.sinc;
  pid_control_V1_B.BusAssignment.state.pose.orientation.x =
    pid_control_V1_B.cosb * pid_control_V1_B.cosc * pid_control_V1_B.sina -
    pid_control_V1_B.cosa * pid_control_V1_B.sinb * pid_control_V1_B.sinc;
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

  /* Gain: '<S43>/ZeroGain' */
  pid_control_V1_B.sina = 0.0 * pid_control_V1_B.SignPreSat_a;

  /* DeadZone: '<S45>/DeadZone' */
  if (pid_control_V1_B.SignPreSat_a > 0.3490658503988659) {
    pid_control_V1_B.SignPreSat_a -= 0.3490658503988659;
  } else if (pid_control_V1_B.SignPreSat_a >= -0.3490658503988659) {
    pid_control_V1_B.SignPreSat_a = 0.0;
  } else {
    pid_control_V1_B.SignPreSat_a -= -0.3490658503988659;
  }

  /* End of DeadZone: '<S45>/DeadZone' */

  /* Gain: '<S50>/Integral Gain' incorporates:
   *  Constant: '<Root>/Constant1'
   *  Sum: '<Root>/Sum4'
   */
  pid_control_V1_B.Switch_n = (0.0 - pid_control_V1_B.x[6]) * -0.01;

  /* Signum: '<S43>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat_a)) {
    /* DataTypeConversion: '<S43>/DataTypeConv1' */
    i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat_a < 0.0) {
      /* DataTypeConversion: '<S43>/DataTypeConv1' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S43>/DataTypeConv1' */
      pid_control_V1_B.chi = (pid_control_V1_B.SignPreSat_a > 0.0);
    }

    /* DataTypeConversion: '<S43>/DataTypeConv1' */
    i = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S43>/SignPreSat' */

  /* Signum: '<S43>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch_n)) {
    /* DataTypeConversion: '<S43>/DataTypeConv2' */
    tmp_1 = 0;
  } else {
    if (pid_control_V1_B.Switch_n < 0.0) {
      /* DataTypeConversion: '<S43>/DataTypeConv2' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S43>/DataTypeConv2' */
      pid_control_V1_B.chi = (pid_control_V1_B.Switch_n > 0.0);
    }

    /* DataTypeConversion: '<S43>/DataTypeConv2' */
    tmp_1 = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S43>/SignPreIntegrator' */

  /* DataTypeConversion: '<S43>/DataTypeConv1' */
  if (i < 0) {
    i = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(i))));
  }

  /* DataTypeConversion: '<S43>/DataTypeConv2' */
  if (tmp_1 < 0) {
    tmp_1 = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(tmp_1))));
  }

  /* Logic: '<S43>/AND3' incorporates:
   *  DataTypeConversion: '<S43>/DataTypeConv1'
   *  DataTypeConversion: '<S43>/DataTypeConv2'
   *  RelationalOperator: '<S43>/Equal1'
   *  RelationalOperator: '<S43>/NotEqual'
   */
  pid_control_V1_B.AND3 = ((pid_control_V1_B.sina !=
    pid_control_V1_B.SignPreSat_a) && (i == tmp_1));
  if (tmp_0) {
    /* Memory: '<S43>/Memory' */
    pid_control_V1_B.Memory_a = pid_control_V1_DW.Memory_PreviousInput_o;
  }

  /* Switch: '<S43>/Switch' */
  if (pid_control_V1_B.Memory_a) {
    /* Gain: '<S50>/Integral Gain' incorporates:
     *  Constant: '<S43>/Constant1'
     *  Switch: '<S43>/Switch'
     */
    pid_control_V1_B.Switch_n = 0.0;
  }

  /* End of Switch: '<S43>/Switch' */

  /* Gain: '<S97>/ZeroGain' */
  pid_control_V1_B.SignPreSat_a = 0.0 * pid_control_V1_B.SignPreSat;

  /* DeadZone: '<S99>/DeadZone' */
  if (pid_control_V1_B.SignPreSat > 20.0) {
    pid_control_V1_B.SignPreSat -= 20.0;
  } else if (pid_control_V1_B.SignPreSat >= 0.0) {
    pid_control_V1_B.SignPreSat = 0.0;
  }

  /* End of DeadZone: '<S99>/DeadZone' */

  /* Gain: '<S104>/Integral Gain' */
  pid_control_V1_B.Switch_k *= 0.01;

  /* Signum: '<S97>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat)) {
    /* DataTypeConversion: '<S97>/DataTypeConv1' */
    i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat < 0.0) {
      /* DataTypeConversion: '<S97>/DataTypeConv1' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S97>/DataTypeConv1' */
      pid_control_V1_B.chi = (pid_control_V1_B.SignPreSat > 0.0);
    }

    /* DataTypeConversion: '<S97>/DataTypeConv1' */
    i = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S97>/SignPreSat' */

  /* Signum: '<S97>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Switch_k)) {
    /* DataTypeConversion: '<S97>/DataTypeConv2' */
    tmp_1 = 0;
  } else {
    if (pid_control_V1_B.Switch_k < 0.0) {
      /* DataTypeConversion: '<S97>/DataTypeConv2' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S97>/DataTypeConv2' */
      pid_control_V1_B.chi = (pid_control_V1_B.Switch_k > 0.0);
    }

    /* DataTypeConversion: '<S97>/DataTypeConv2' */
    tmp_1 = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S97>/SignPreIntegrator' */

  /* DataTypeConversion: '<S97>/DataTypeConv1' */
  if (i < 0) {
    i = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(i))));
  }

  /* DataTypeConversion: '<S97>/DataTypeConv2' */
  if (tmp_1 < 0) {
    tmp_1 = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(tmp_1))));
  }

  /* Logic: '<S97>/AND3' incorporates:
   *  DataTypeConversion: '<S97>/DataTypeConv1'
   *  DataTypeConversion: '<S97>/DataTypeConv2'
   *  RelationalOperator: '<S97>/Equal1'
   *  RelationalOperator: '<S97>/NotEqual'
   */
  pid_control_V1_B.AND3_e = ((pid_control_V1_B.SignPreSat_a !=
    pid_control_V1_B.SignPreSat) && (i == tmp_1));
  if (tmp_0) {
    /* Memory: '<S97>/Memory' */
    pid_control_V1_B.Memory_n = pid_control_V1_DW.Memory_PreviousInput_m;
  }

  /* Switch: '<S97>/Switch' */
  if (pid_control_V1_B.Memory_n) {
    /* Gain: '<S104>/Integral Gain' incorporates:
     *  Constant: '<S97>/Constant1'
     *  Switch: '<S97>/Switch'
     */
    pid_control_V1_B.Switch_k = 0.0;
  }

  /* End of Switch: '<S97>/Switch' */

  /* Sum: '<S151>/SumI4' incorporates:
   *  Gain: '<S156>/Integral Gain'
   *  Sum: '<S151>/SumI2'
   */
  pid_control_V1_B.SumI4 = (pid_control_V1_B.Saturation_f -
    pid_control_V1_B.Sum_hl) + -0.5 * pid_control_V1_B.Sum1_g;

  /* Gain: '<S208>/Integral Gain' */
  pid_control_V1_B.IntegralGain = -0.05 * pid_control_V1_B.Sum5;

  /* Gain: '<S268>/Filter Coefficient' incorporates:
   *  Constant: '<Root>/Constant3'
   *  Gain: '<S258>/Derivative Gain'
   *  Integrator: '<S260>/Filter'
   *  Sum: '<Root>/Sum3'
   *  Sum: '<S260>/SumD'
   */
  pid_control_V1_B.FilterCoefficient_cv = ((28.0 - pid_control_V1_B.x[0]) *
    0.005 - pid_control_V1_X.Filter_CSTATE_l) * 100.0;

  /* Sum: '<S274>/Sum' incorporates:
   *  Constant: '<Root>/Constant3'
   *  Gain: '<S270>/Proportional Gain'
   *  Integrator: '<S265>/Integrator'
   *  Sum: '<Root>/Sum3'
   */
  pid_control_V1_B.Sum1_g = ((28.0 - pid_control_V1_B.x[0]) * 0.05 +
    pid_control_V1_X.Integrator_CSTATE_f) +
    pid_control_V1_B.FilterCoefficient_cv;

  /* DeadZone: '<S257>/DeadZone' */
  if (pid_control_V1_B.Sum1_g > 1.0) {
    pid_control_V1_B.SignPreSat = pid_control_V1_B.Sum1_g - 1.0;
  } else if (pid_control_V1_B.Sum1_g >= 0.0) {
    pid_control_V1_B.SignPreSat = 0.0;
  } else {
    pid_control_V1_B.SignPreSat = pid_control_V1_B.Sum1_g;
  }

  /* End of DeadZone: '<S257>/DeadZone' */

  /* Gain: '<S262>/Integral Gain' incorporates:
   *  Constant: '<Root>/Constant3'
   *  Sum: '<Root>/Sum3'
   */
  pid_control_V1_B.Sum_hl = (28.0 - pid_control_V1_B.x[0]) * 0.01;

  /* Signum: '<S255>/SignPreSat' */
  if (rtIsNaN(pid_control_V1_B.SignPreSat)) {
    /* DataTypeConversion: '<S255>/DataTypeConv1' */
    i = 0;
  } else {
    if (pid_control_V1_B.SignPreSat < 0.0) {
      /* DataTypeConversion: '<S255>/DataTypeConv1' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S255>/DataTypeConv1' */
      pid_control_V1_B.chi = (pid_control_V1_B.SignPreSat > 0.0);
    }

    /* DataTypeConversion: '<S255>/DataTypeConv1' */
    i = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S255>/SignPreSat' */

  /* Signum: '<S255>/SignPreIntegrator' */
  if (rtIsNaN(pid_control_V1_B.Sum_hl)) {
    /* DataTypeConversion: '<S255>/DataTypeConv2' */
    tmp_1 = 0;
  } else {
    if (pid_control_V1_B.Sum_hl < 0.0) {
      /* DataTypeConversion: '<S255>/DataTypeConv2' */
      pid_control_V1_B.chi = -1.0;
    } else {
      /* DataTypeConversion: '<S255>/DataTypeConv2' */
      pid_control_V1_B.chi = (pid_control_V1_B.Sum_hl > 0.0);
    }

    /* DataTypeConversion: '<S255>/DataTypeConv2' */
    tmp_1 = static_cast<int32_T>(fmod(pid_control_V1_B.chi, 256.0));
  }

  /* End of Signum: '<S255>/SignPreIntegrator' */

  /* DataTypeConversion: '<S255>/DataTypeConv1' */
  if (i < 0) {
    i = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(i))));
  }

  /* DataTypeConversion: '<S255>/DataTypeConv2' */
  if (tmp_1 < 0) {
    tmp_1 = static_cast<int8_T>(-static_cast<int8_T>(static_cast<uint8_T>(-
      static_cast<real_T>(tmp_1))));
  }

  /* Logic: '<S255>/AND3' incorporates:
   *  DataTypeConversion: '<S255>/DataTypeConv1'
   *  DataTypeConversion: '<S255>/DataTypeConv2'
   *  Gain: '<S255>/ZeroGain'
   *  RelationalOperator: '<S255>/Equal1'
   *  RelationalOperator: '<S255>/NotEqual'
   */
  pid_control_V1_B.AND3_c = ((0.0 * pid_control_V1_B.Sum1_g !=
    pid_control_V1_B.SignPreSat) && (i == tmp_1));
  if (tmp_0) {
    /* Memory: '<S255>/Memory' */
    pid_control_V1_B.Memory_h = pid_control_V1_DW.Memory_PreviousInput_a;
  }

  /* Switch: '<S255>/Switch' */
  if (pid_control_V1_B.Memory_h) {
    /* Switch: '<S255>/Switch' incorporates:
     *  Constant: '<S255>/Constant1'
     */
    pid_control_V1_B.Switch_j = 0.0;
  } else {
    /* Switch: '<S255>/Switch' */
    pid_control_V1_B.Switch_j = pid_control_V1_B.Sum_hl;
  }

  /* End of Switch: '<S255>/Switch' */

  /* Saturate: '<S272>/Saturation' */
  if (pid_control_V1_B.Sum1_g > 1.0) {
    /* Saturate: '<S272>/Saturation' */
    pid_control_V1_B.Saturation_o = 1.0;
  } else if (pid_control_V1_B.Sum1_g < 0.0) {
    /* Saturate: '<S272>/Saturation' */
    pid_control_V1_B.Saturation_o = 0.0;
  } else {
    /* Saturate: '<S272>/Saturation' */
    pid_control_V1_B.Saturation_o = pid_control_V1_B.Sum1_g;
  }

  /* End of Saturate: '<S272>/Saturation' */
  if (tmp_0) {
    /* Memory: '<S16>/Memory' */
    pid_control_V1_B.Memory[0] = pid_control_V1_DW.Memory_PreviousInput[0];

    /* Memory: '<S16>/Memory1' */
    pid_control_V1_B.Memory1[0] = pid_control_V1_DW.Memory1_PreviousInput[0];

    /* Memory: '<S16>/Memory' */
    pid_control_V1_B.Memory[1] = pid_control_V1_DW.Memory_PreviousInput[1];

    /* Memory: '<S16>/Memory1' */
    pid_control_V1_B.Memory1[1] = pid_control_V1_DW.Memory1_PreviousInput[1];

    /* Memory: '<S16>/Memory' */
    pid_control_V1_B.Memory[2] = pid_control_V1_DW.Memory_PreviousInput[2];

    /* Memory: '<S16>/Memory1' */
    pid_control_V1_B.Memory1[2] = pid_control_V1_DW.Memory1_PreviousInput[2];
  }

  /* SignalConversion generated from: '<S290>/ SFunction ' incorporates:
   *  MATLAB Function: '<S16>/MATLAB Function'
   */
  pid_control_V1_B.TmpSignalConversionAtSFunct[0] =
    pid_control_V1_B.Saturation_k;
  pid_control_V1_B.TmpSignalConversionAtSFunct[1] =
    pid_control_V1_B.Saturation_f;
  pid_control_V1_B.TmpSignalConversionAtSFunct[2] =
    pid_control_V1_B.Saturation_m;
  pid_control_V1_B.TmpSignalConversionAtSFunct[3] =
    pid_control_V1_B.Saturation_o;
  pid_control_V1_B.TmpSignalConversionAtSFunct[4] =
    pid_control_V1_B.Saturation_o;

  /* MATLAB Function: '<S16>/MATLAB Function' incorporates:
   *  Memory: '<S16>/Memory'
   */
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[1] <= 0.3490658503988659) {
    pid_control_V1_B.u2 = pid_control_V1_B.TmpSignalConversionAtSFunct[1];
  } else {
    pid_control_V1_B.u2 = 0.3490658503988659;
  }

  if (!(pid_control_V1_B.u2 >= -0.3490658503988659)) {
    pid_control_V1_B.u2 = -0.3490658503988659;
  }

  tmp_2 = _mm_add_pd(_mm_loadu_pd(&pid_control_V1_B.x[0]), _mm_loadu_pd
                     (&pid_control_V1_B.Memory[0]));
  _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);

  /* MATLAB Function: '<S16>/MATLAB Function' incorporates:
   *  Memory: '<S16>/Memory'
   *  Memory: '<S16>/Memory1'
   */
  pid_control_V1_B.SignPreSat = pid_control_V1_B.x[2] + pid_control_V1_B.Memory
    [2];
  pid_control_V1_B.chi = sqrt((pid_control_V1_B.dv[0] * pid_control_V1_B.dv[0] +
    pid_control_V1_B.dv[1] * pid_control_V1_B.dv[1]) +
    pid_control_V1_B.SignPreSat * pid_control_V1_B.SignPreSat);
  if (pid_control_V1_B.chi == 0.0) {
    pid_control_V1_B.chi = 0.001;
  }

  pid_control_V1_B.SignPreSat = rt_atan2d_snf(pid_control_V1_B.SignPreSat,
    pid_control_V1_B.dv[0]);
  pid_control_V1_B.SignPreSat_a = asin(pid_control_V1_B.dv[1] /
    pid_control_V1_B.chi);
  if ((-pid_control_V1_B.x[11] - 0.363 <= 0.001) || rtIsNaN(-pid_control_V1_B.x
       [11] - 0.363)) {
    pid_control_V1_B.Sum_hl = 0.001;
  } else {
    pid_control_V1_B.Sum_hl = -pid_control_V1_B.x[11] - 0.363;
  }

  if ((-pid_control_V1_B.x[11] + 2.5 <= 0.001) || rtIsNaN(-pid_control_V1_B.x[11]
       + 2.5)) {
    pid_control_V1_B.Sum1_g = 0.001;
  } else {
    pid_control_V1_B.Sum1_g = -pid_control_V1_B.x[11] + 2.5;
  }

  pid_control_V1_B.Q = pid_control_V1_B.chi * pid_control_V1_B.chi * 0.6125;
  pid_control_V1_B.wbe_b[0] = pid_control_V1_B.x[3];
  pid_control_V1_B.wbe_b[1] = pid_control_V1_B.x[4];
  pid_control_V1_B.wbe_b[2] = pid_control_V1_B.x[5];
  pid_control_V1_B.sina = ((pid_control_V1_B.SignPreSat - -0.065449846949787352)
    + 0.017453292519943295) * 4.9604094530365153;
  pid_control_V1_B.sinb = ((pid_control_V1_B.SignPreSat - -0.074176493209759012)
    + 0.0087266462599716477) * 4.8387748917360032;
  pid_control_V1_B.Cn = pid_control_V1_B.Sum_hl / 5.02;
  pid_control_V1_B.sinc = (rt_powd_snf(pid_control_V1_B.Cn, 0.787) * 288.0 * exp
    (rt_powd_snf(pid_control_V1_B.Cn, 0.327) * -9.14) * 0.97986308862072491 /
    5.9129476540958859 + 1.0) * pid_control_V1_B.sina;
  pid_control_V1_B.Sum1_g /= 2.74;
  pid_control_V1_B.cosa = (rt_powd_snf(pid_control_V1_B.Sum1_g, 0.787) * 288.0 *
    exp(rt_powd_snf(pid_control_V1_B.Sum1_g, 0.327) * -9.14) *
    0.95628590200128227 / 5.35300902982722 + 1.0) * pid_control_V1_B.sinb;
  pid_control_V1_B.cosb = (1.0 - exp(rt_powd_snf(pid_control_V1_B.Cn, 0.686) *
    -10.1)) * (pid_control_V1_B.sina * pid_control_V1_B.sina /
               21.205750411731103);
  pid_control_V1_B.cosc = (1.0 - exp(rt_powd_snf(pid_control_V1_B.Sum1_g, 0.686)
    * -10.1)) * (pid_control_V1_B.sinb * pid_control_V1_B.sinb /
                 18.943803701146454);
  pid_control_V1_B.Dtot = ((pid_control_V1_B.u2 * pid_control_V1_B.u2 * -1.08E-5
    + 0.000715 * pid_control_V1_B.u2) * 1.128 + ((pid_control_V1_B.cosb * 3.334
    + 0.1020204) + pid_control_V1_B.cosc * 1.128)) * pid_control_V1_B.Q;
  pid_control_V1_B.Ltot = (pid_control_V1_B.sinc * 3.334 + pid_control_V1_B.cosa
    * 1.128) * pid_control_V1_B.Q;
  pid_control_V1_B.CQ = -0.019 * pid_control_V1_B.SignPreSat_a * 180.0 /
    3.1415926535897931;
  pid_control_V1_B.Sum1_g = sin(pid_control_V1_B.SignPreSat);
  pid_control_V1_B.Sum_hl = cos(pid_control_V1_B.SignPreSat);
  pid_control_V1_B.FA_b_tmp[0] = pid_control_V1_B.Sum_hl;
  pid_control_V1_B.FA_b_tmp[3] = 0.0;
  pid_control_V1_B.FA_b_tmp[6] = -pid_control_V1_B.Sum1_g;
  pid_control_V1_B.FA_b_tmp[2] = pid_control_V1_B.Sum1_g;
  pid_control_V1_B.FA_b_tmp[5] = 0.0;
  pid_control_V1_B.FA_b_tmp[8] = pid_control_V1_B.Sum_hl;
  pid_control_V1_B.FA_b[0] = -pid_control_V1_B.Dtot;
  pid_control_V1_B.FA_b[1] = pid_control_V1_B.CQ * pid_control_V1_B.Q * 3.334;
  pid_control_V1_B.FA_b[2] = -pid_control_V1_B.Ltot;
  pid_control_V1_B.FA_b_tmp[1] = 0.0;
  pid_control_V1_B.FA_b_tmp[4] = 1.0;
  pid_control_V1_B.FA_b_tmp[7] = 0.0;
  pid_control_V1_B.Sum1_g = 0.0;
  pid_control_V1_B.Sum_hl = 0.0;
  pid_control_V1_B.Sum5 = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.FA_b[i])), _mm_set_pd(pid_control_V1_B.Sum_hl,
      pid_control_V1_B.Sum1_g));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
    pid_control_V1_B.Sum1_g = pid_control_V1_B.dv[0];
    pid_control_V1_B.Sum_hl = pid_control_V1_B.dv[1];
    pid_control_V1_B.Sum5 += pid_control_V1_B.FA_b_tmp[3 * i + 2] *
      pid_control_V1_B.FA_b[i];
  }

  if (pid_control_V1_B.TmpSignalConversionAtSFunct[0] <= 0.26179938779914941) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[0];
  } else {
    pid_control_V1_B.Vd1 = 0.26179938779914941;
  }

  if (!(pid_control_V1_B.Vd1 >= -0.3490658503988659)) {
    pid_control_V1_B.Vd1 = -0.3490658503988659;
  }

  pid_control_V1_B.FE1_b_idx_1 = 2.0 * pid_control_V1_B.chi;
  pid_control_V1_B.Cl = ((pid_control_V1_B.Memory1[0] + pid_control_V1_B.x[3]) *
    5.02 / pid_control_V1_B.FE1_b_idx_1 * -2.0 + -0.0286 *
    pid_control_V1_B.SignPreSat_a) + -0.5 * pid_control_V1_B.Vd1;
  pid_control_V1_B.u2 = ((pid_control_V1_B.Memory1[1] + pid_control_V1_B.x[4]) *
    0.646 / pid_control_V1_B.FE1_b_idx_1 * -5.0 + (exp(pid_control_V1_B.Cn *
    -4.0) * -0.05 + -1.14 * pid_control_V1_B.SignPreSat)) + -3.0 *
    pid_control_V1_B.u2;
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[2] <= 0.26179938779914941) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[2];
  } else {
    pid_control_V1_B.Vd1 = 0.26179938779914941;
  }

  if (!(pid_control_V1_B.Vd1 >= -0.26179938779914941)) {
    pid_control_V1_B.Vd1 = -0.26179938779914941;
  }

  pid_control_V1_B.Cn = ((pid_control_V1_B.Memory1[2] + pid_control_V1_B.x[5]) *
    5.02 / pid_control_V1_B.FE1_b_idx_1 * -1.5 + -0.1146 *
    pid_control_V1_B.SignPreSat_a) + -0.3 * pid_control_V1_B.Vd1;
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[3] <= 1.0) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[3];
  } else {
    pid_control_V1_B.Vd1 = 1.0;
  }

  if (!(pid_control_V1_B.Vd1 >= 0.0)) {
    pid_control_V1_B.Vd1 = 0.0;
  }

  pid_control_V1_B.Vd1 = (37.42 - pid_control_V1_B.chi) * pid_control_V1_B.Vd1 +
    pid_control_V1_B.chi;
  pid_control_V1_B.Tp1 = 0.11056515539994617 * pid_control_V1_B.Vd1 *
    (pid_control_V1_B.Vd1 - pid_control_V1_B.chi);
  if (pid_control_V1_B.TmpSignalConversionAtSFunct[4] <= 1.0) {
    pid_control_V1_B.Vd1 = pid_control_V1_B.TmpSignalConversionAtSFunct[4];
  } else {
    pid_control_V1_B.Vd1 = 1.0;
  }

  if (!(pid_control_V1_B.Vd1 >= 0.0)) {
    pid_control_V1_B.Vd1 = 0.0;
  }

  pid_control_V1_B.Vd1 = (37.42 - pid_control_V1_B.chi) * pid_control_V1_B.Vd1 +
    pid_control_V1_B.chi;
  pid_control_V1_B.Tp2 = 0.11056515539994617 * pid_control_V1_B.Vd1 *
    (pid_control_V1_B.Vd1 - pid_control_V1_B.chi);
  pid_control_V1_B.Vd1 = pid_control_V1_B.Tp1 * 0.99619469809174555;
  pid_control_V1_B.FE1_b_idx_1 = 0.0;
  pid_control_V1_B.Tp1 *= 0.087155742747658166;
  pid_control_V1_B.FE2_b_idx_0 = pid_control_V1_B.Tp2 * 0.99619469809174555;
  pid_control_V1_B.FE2_b_idx_2 = pid_control_V1_B.Tp2 * 0.087155742747658166;
  pid_control_V1_B.Tp2 = -9.81 * sin(pid_control_V1_B.x[7]) * 112.0;
  _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_mul_pd(_mm_mul_pd(_mm_mul_pd
    (_mm_set1_pd(9.81), _mm_set1_pd(cos(pid_control_V1_B.x[7]))), _mm_set_pd(cos
    (pid_control_V1_B.x[6]), sin(pid_control_V1_B.x[6]))), _mm_set1_pd(112.0)));
  pid_control_V1_B.Fg_b_idx_1 = pid_control_V1_B.dv[0];
  pid_control_V1_B.Fg_b_idx_2 = pid_control_V1_B.dv[1];
  pid_control_V1_B.Mcg_b_idx_2 = 5.02 * pid_control_V1_B.Q * 3.334;
  pid_control_V1_B.Mcg_b_idx_0 = (0.6 * pid_control_V1_B.Tp1 + -0.6 *
    pid_control_V1_B.FE2_b_idx_2) + pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cl;
  pid_control_V1_B.Q = 0.646 * pid_control_V1_B.Q * 3.334 * pid_control_V1_B.u2
    + ((-0.285 * pid_control_V1_B.Vd1 - 0.519 * pid_control_V1_B.Tp1) + (-0.285 *
        pid_control_V1_B.FE2_b_idx_0 - 0.519 * pid_control_V1_B.FE2_b_idx_2));
  pid_control_V1_B.Mcg_b_idx_2 = ((0.0 - 0.6 * pid_control_V1_B.Vd1) + (0.0 -
    -0.6 * pid_control_V1_B.FE2_b_idx_0)) + pid_control_V1_B.Mcg_b_idx_2 *
    pid_control_V1_B.Cn;
  pid_control_V1_B.c_phi = cos(pid_control_V1_B.x[6]);
  pid_control_V1_B.s_phi = sin(pid_control_V1_B.x[6]);
  pid_control_V1_B.c_the = cos(pid_control_V1_B.x[7]);
  pid_control_V1_B.s_the = sin(pid_control_V1_B.x[7]);
  pid_control_V1_B.c_psi = cos(pid_control_V1_B.x[8]);
  pid_control_V1_B.s_psi = sin(pid_control_V1_B.x[8]);
  pid_control_V1_B.FA_b_tmp[0] = pid_control_V1_B.c_the * pid_control_V1_B.c_psi;
  pid_control_V1_B.c_the_tmp = pid_control_V1_B.s_phi * pid_control_V1_B.s_the;
  pid_control_V1_B.FA_b_tmp[3] = pid_control_V1_B.c_the_tmp *
    pid_control_V1_B.c_psi - pid_control_V1_B.c_phi * pid_control_V1_B.s_psi;
  pid_control_V1_B.c_the_tmp_g = pid_control_V1_B.c_phi * pid_control_V1_B.s_the;
  pid_control_V1_B.FA_b_tmp[6] = pid_control_V1_B.c_the_tmp_g *
    pid_control_V1_B.c_psi + pid_control_V1_B.s_phi * pid_control_V1_B.s_psi;
  pid_control_V1_B.FA_b_tmp[1] = pid_control_V1_B.c_the * pid_control_V1_B.s_psi;
  pid_control_V1_B.FA_b_tmp[4] = pid_control_V1_B.c_the_tmp *
    pid_control_V1_B.s_psi + pid_control_V1_B.c_phi * pid_control_V1_B.c_psi;
  pid_control_V1_B.FA_b_tmp[7] = pid_control_V1_B.c_the_tmp_g *
    pid_control_V1_B.s_psi - pid_control_V1_B.s_phi * pid_control_V1_B.c_psi;
  pid_control_V1_B.FA_b_tmp[2] = -pid_control_V1_B.s_the;
  pid_control_V1_B.FA_b_tmp[5] = pid_control_V1_B.s_phi * pid_control_V1_B.c_the;
  pid_control_V1_B.FA_b_tmp[8] = pid_control_V1_B.c_phi * pid_control_V1_B.c_the;
  pid_control_V1_B.wbe_b_m[0] = pid_control_V1_B.x[0];
  pid_control_V1_B.wbe_b_m[1] = pid_control_V1_B.x[1];
  pid_control_V1_B.wbe_b_m[2] = pid_control_V1_B.x[2];
  pid_control_V1_B.FE2_b_idx_0 += pid_control_V1_B.Vd1;
  pid_control_V1_B.c_phi = pid_control_V1_B.FE2_b_idx_0;
  pid_control_V1_B.F_b[0] = (pid_control_V1_B.Tp2 + pid_control_V1_B.FE2_b_idx_0)
    + pid_control_V1_B.Sum1_g;
  pid_control_V1_B.Vd1 = 0.0;
  pid_control_V1_B.F_b[1] = pid_control_V1_B.dv[0] + pid_control_V1_B.Sum_hl;
  pid_control_V1_B.FE2_b_idx_0 = pid_control_V1_B.Tp1 +
    pid_control_V1_B.FE2_b_idx_2;
  pid_control_V1_B.F_b[2] = (pid_control_V1_B.dv[1] +
    pid_control_V1_B.FE2_b_idx_0) + pid_control_V1_B.Sum5;
  pid_control_V1_B.Tp1 = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.wbe_b_m[i])), _mm_set_pd
                       (pid_control_V1_B.FE1_b_idx_1, pid_control_V1_B.Vd1));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
    pid_control_V1_B.Vd1 = pid_control_V1_B.dv[0];
    pid_control_V1_B.FE1_b_idx_1 = pid_control_V1_B.dv[1];
    pid_control_V1_B.Tp1 += pid_control_V1_B.FA_b_tmp[3 * i + 2] *
      pid_control_V1_B.wbe_b_m[i];
  }

  tmp_2 = _mm_sub_pd(_mm_mul_pd(_mm_set_pd(pid_control_V1_B.x[0],
    pid_control_V1_B.x[2]), _mm_loadu_pd(&pid_control_V1_B.x[4])), _mm_mul_pd
                     (_mm_loadu_pd(&pid_control_V1_B.x[1]), _mm_set_pd
                      (pid_control_V1_B.x[3], pid_control_V1_B.x[5])));
  _mm_storeu_pd(&pid_control_V1_B.wbe_b_m[0], tmp_2);
  pid_control_V1_B.wbe_b_m[2] = pid_control_V1_B.x[1] * pid_control_V1_B.x[3] -
    pid_control_V1_B.x[0] * pid_control_V1_B.x[4];
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
  pid_control_V1_B.FE2_b_idx_2 = 0.0;
  pid_control_V1_B.s_phi = 0.0;
  pid_control_V1_B.c_the = 0.0;
  for (i = 0; i < 3; i++) {
    tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&pid_control_V1_B.FA_b_tmp[3 * i]),
      _mm_set1_pd(pid_control_V1_B.wbe_b[i])), _mm_set_pd(pid_control_V1_B.s_phi,
      pid_control_V1_B.FE2_b_idx_2));
    _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
    pid_control_V1_B.FE2_b_idx_2 = pid_control_V1_B.dv[0];
    pid_control_V1_B.s_phi = pid_control_V1_B.dv[1];
    _mm_storeu_pd(&pid_control_V1_B.dv[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (0.0089285714285714281, pid_control_V1_B.FA_b_tmp[3 * i + 2]), _mm_set_pd
      (pid_control_V1_B.F_b[i], pid_control_V1_B.wbe_b[i])), _mm_mul_pd
      (_mm_set_pd(pid_control_V1_B.wbe_b_m[i], pid_control_V1_B.c_the),
       _mm_set_pd(-1.0, 1.0))));
    pid_control_V1_B.c_the = pid_control_V1_B.dv[0];
    pid_control_V1_B.XDOT[i] = pid_control_V1_B.dv[1];
  }

  pid_control_V1_B.XDOT[3] = ((114.39 * pid_control_V1_B.Mcg_b_idx_0 + 8.97 *
    pid_control_V1_B.Mcg_b_idx_2) - (-615.2523000000001 * pid_control_V1_B.x[3]
    + 3384.0440999999996 * pid_control_V1_B.x[5]) * pid_control_V1_B.x[4]) /
    4461.966;
  pid_control_V1_B.XDOT[4] = ((pid_control_V1_B.Q - -74.68 * pid_control_V1_B.x
    [3] * pid_control_V1_B.x[5]) - (pid_control_V1_B.x[3] * pid_control_V1_B.x[3]
    - pid_control_V1_B.x[5] * pid_control_V1_B.x[5]) * 8.97) / 85.51;
  pid_control_V1_B.XDOT[5] = ((-615.2523000000001 * pid_control_V1_B.x[5] +
    -1738.2571000000003 * pid_control_V1_B.x[3]) * pid_control_V1_B.x[4] + (8.97
    * pid_control_V1_B.Mcg_b_idx_0 + 39.71 * pid_control_V1_B.Mcg_b_idx_2)) /
    4461.966;
  pid_control_V1_B.XDOT[9] = pid_control_V1_B.Vd1;
  pid_control_V1_B.XDOT[10] = pid_control_V1_B.FE1_b_idx_1;
  pid_control_V1_B.XDOT[11] = pid_control_V1_B.Tp1;
  pid_control_V1_B.XDOT[12] = pid_control_V1_B.Ltot / pid_control_V1_B.Dtot;
  pid_control_V1_B.XDOT[19] = pid_control_V1_B.CQ;
  pid_control_V1_B.XDOT[20] = pid_control_V1_B.Cl;
  pid_control_V1_B.XDOT[21] = pid_control_V1_B.u2;
  pid_control_V1_B.XDOT[22] = pid_control_V1_B.Cn;
  pid_control_V1_B.XDOT[23] = pid_control_V1_B.SignPreSat;
  pid_control_V1_B.XDOT[24] = pid_control_V1_B.SignPreSat_a;
  pid_control_V1_B.XDOT[25] = pid_control_V1_B.sina;
  pid_control_V1_B.XDOT[26] = pid_control_V1_B.sinb;
  pid_control_V1_B.XDOT[27] = pid_control_V1_B.sinc;
  pid_control_V1_B.XDOT[28] = pid_control_V1_B.cosa;
  pid_control_V1_B.XDOT[29] = pid_control_V1_B.cosb;
  pid_control_V1_B.XDOT[30] = pid_control_V1_B.cosc;
  pid_control_V1_B.XDOT[6] = pid_control_V1_B.FE2_b_idx_2;
  pid_control_V1_B.XDOT[13] = pid_control_V1_B.F_b[0];
  pid_control_V1_B.XDOT[16] = pid_control_V1_B.Mcg_b_idx_0;
  pid_control_V1_B.XDOT[31] = pid_control_V1_B.Tp2;
  pid_control_V1_B.XDOT[34] = pid_control_V1_B.c_phi;
  pid_control_V1_B.XDOT[37] = pid_control_V1_B.Sum1_g;
  pid_control_V1_B.XDOT[7] = pid_control_V1_B.s_phi;
  pid_control_V1_B.XDOT[14] = pid_control_V1_B.F_b[1];
  pid_control_V1_B.XDOT[17] = pid_control_V1_B.Q;
  pid_control_V1_B.XDOT[32] = pid_control_V1_B.Fg_b_idx_1;
  pid_control_V1_B.XDOT[35] = 0.0;
  pid_control_V1_B.XDOT[38] = pid_control_V1_B.Sum_hl;
  pid_control_V1_B.XDOT[8] = pid_control_V1_B.c_the;
  pid_control_V1_B.XDOT[15] = pid_control_V1_B.F_b[2];
  pid_control_V1_B.XDOT[18] = pid_control_V1_B.Mcg_b_idx_2;
  pid_control_V1_B.XDOT[33] = pid_control_V1_B.Fg_b_idx_2;
  pid_control_V1_B.XDOT[36] = pid_control_V1_B.FE2_b_idx_0;
  pid_control_V1_B.XDOT[39] = pid_control_V1_B.Sum5;

  /* Math: '<S16>/Square2' */
  pid_control_V1_B.Sum1_g = pid_control_V1_B.x[1] * pid_control_V1_B.x[1];

  /* Product: '<S16>/Product2' incorporates:
   *  Math: '<S16>/Square'
   *  Math: '<S16>/Square1'
   *  Sqrt: '<S16>/Sqrt'
   *  Sum: '<S16>/Sum2'
   */
  pid_control_V1_B.Power = sqrt((pid_control_V1_B.x[0] * pid_control_V1_B.x[0] +
    pid_control_V1_B.Sum1_g) + pid_control_V1_B.x[2] * pid_control_V1_B.x[2]) *
    pid_control_V1_B.XDOT[34];

  /* Gain: '<S16>/Gain3' */
  pid_control_V1_B.Gain3 = 0.001 * pid_control_V1_B.Power;
  if (tmp_0) {
  }

  /* Gain: '<S16>/Gain1' incorporates:
   *  Integrator: '<S16>/Integrator1'
   */
  pid_control_V1_B.EnergykWh = 2.7777777777777776E-7 *
    pid_control_V1_X.Integrator1_CSTATE;
  if (tmp_0) {
    /* SignalConversion generated from: '<S16>/To Workspace1' */
    pid_control_V1_B.TmpSignalConversionAtSFunct[0] =
      pid_control_V1_B.Saturation_k;
    pid_control_V1_B.TmpSignalConversionAtSFunct[1] =
      pid_control_V1_B.Saturation_f;
    pid_control_V1_B.TmpSignalConversionAtSFunct[2] =
      pid_control_V1_B.Saturation_m;
    pid_control_V1_B.TmpSignalConversionAtSFunct[3] =
      pid_control_V1_B.Saturation_o;
    pid_control_V1_B.TmpSignalConversionAtSFunct[4] =
      pid_control_V1_B.Saturation_o;
  }

  /* Product: '<S16>/Divide' incorporates:
   *  Constant: '<S16>/thrust efficiency Cp?'
   */
  pid_control_V1_B.powerdemand = pid_control_V1_B.Gain3 / 0.248;
  if (tmp_0) {
  }

  /* Product: '<S16>/Divide1' */
  pid_control_V1_B.loadtorque = pid_control_V1_B.powerdemand /
    pid_control_V1_ConstB.motorspeed;
  if (tmp_0) {
    /* Gain: '<S288>/Output' incorporates:
     *  RandomNumber: '<S288>/White Noise'
     */
    pid_control_V1_B.Output = 10.0 * pid_control_V1_DW.NextOutput;
  }

  /* UnitConversion: '<S298>/Unit Conversion' incorporates:
   *  Gain: '<S16>/Gain4'
   */
  /* Unit Conversion - from: m to: ft
     Expression: output = (3.28084*input) + (0) */
  pid_control_V1_B.SignPreSat = 3.280839895013123 * -pid_control_V1_B.x[11];

  /* Saturate: '<S331>/Limit Function 10ft to 1000ft' */
  if (pid_control_V1_B.SignPreSat > 1000.0) {
    pid_control_V1_B.Sum_hl = 1000.0;
  } else if (pid_control_V1_B.SignPreSat < 10.0) {
    pid_control_V1_B.Sum_hl = 10.0;
  } else {
    pid_control_V1_B.Sum_hl = pid_control_V1_B.SignPreSat;
  }

  /* End of Saturate: '<S331>/Limit Function 10ft to 1000ft' */

  /* Gain: '<S303>/Lw' */
  pid_control_V1_B.Sum5 = pid_control_V1_ConstB.UnitConversion_c;

  /* Interpolation_n-D: '<S313>/Medium//High Altitude Intensity' incorporates:
   *  PreLookup: '<S313>/PreLook-Up Index Search  (altitude)'
   */
  pid_control_V1_B.bpIndex[0] = plook_bincpa(pid_control_V1_B.SignPreSat,
    pid_control_V1_ConstP.PreLookUpIndexSearchaltitude_Br, 11U,
    &pid_control_V1_B.Sum1_g, &pid_control_V1_DW.PreLookUpIndexSearchaltitude_DW);
  pid_control_V1_B.frac[0] = pid_control_V1_B.Sum1_g;
  pid_control_V1_B.frac[1] = pid_control_V1_ConstB.PreLookUpIndexSearchprobofe;
  pid_control_V1_B.bpIndex[1] =
    pid_control_V1_ConstB.PreLookUpIndexSearchprobo_g;
  pid_control_V1_B.Sum1_g = intrp2d_la_pw(pid_control_V1_B.bpIndex,
    pid_control_V1_B.frac, pid_control_V1_ConstP.MediumHighAltitudeIntensity_Tab,
    12U, pid_control_V1_ConstP.MediumHighAltitudeIntensity_max);

  /* UnitConversion: '<S304>/Unit Conversion' incorporates:
   *  MATLAB Function: '<S16>/MATLAB Function'
   */
  /* Unit Conversion - from: m/s to: ft/s
     Expression: output = (3.28084*input) + (0) */
  pid_control_V1_B.SignPreSat_a = 3.280839895013123 * pid_control_V1_B.chi;

  /* Outputs for Enabled SubSystem: '<S296>/Hpgw' incorporates:
   *  EnablePort: '<S307>/Enable'
   */
  if (tmp_0) {
    /* Product: '<S306>/Divide' incorporates:
     *  Product: '<S306>/Product'
     *  RandomNumber: '<S306>/White Noise'
     */
    tmp_2 = _mm_mul_pd(_mm_loadu_pd(&pid_control_V1_ConstB.Divide[0]),
                       _mm_loadu_pd(&pid_control_V1_DW.NextOutput_j[0]));

    /* RandomNumber: '<S306>/White Noise' incorporates:
     *  Product: '<S306>/Product'
     */
    _mm_storeu_pd(&pid_control_V1_B.Product[0], tmp_2);

    /* Product: '<S306>/Divide' incorporates:
     *  Product: '<S306>/Product'
     *  RandomNumber: '<S306>/White Noise'
     */
    tmp_2 = _mm_mul_pd(_mm_loadu_pd(&pid_control_V1_ConstB.Divide[2]),
                       _mm_loadu_pd(&pid_control_V1_DW.NextOutput_j[2]));

    /* RandomNumber: '<S306>/White Noise' incorporates:
     *  Product: '<S306>/Product'
     */
    _mm_storeu_pd(&pid_control_V1_B.Product[2], tmp_2);
    if (rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
        (!pid_control_V1_DW.Hpgw_MODE)) {
      (void) memset(&(pid_control_V1_XDis.pgw_p_CSTATE), 0,
                    2*sizeof(boolean_T));

      /* InitializeConditions for Integrator: '<S307>/pgw_p' */
      pid_control_V1_X.pgw_p_CSTATE[0] = 0.0;
      pid_control_V1_X.pgw_p_CSTATE[1] = 0.0;
      pid_control_V1_DW.Hpgw_MODE = true;
    }
  }

  if (pid_control_V1_DW.Hpgw_MODE) {
    /* Fcn: '<S307>/sqrt(0.8//V)' */
    pid_control_V1_B.chi = sqrt(0.8 / pid_control_V1_B.SignPreSat_a);

    /* Product: '<S307>/w3' */
    pid_control_V1_B.sina = pid_control_V1_B.SignPreSat_a *
      pid_control_V1_ConstB.w4;

    /* Product: '<S307>/w' incorporates:
     *  Fcn: '<S307>/sqrt(0.8//V)'
     *  Gain: '<S303>/Lw'
     *  Integrator: '<S307>/pgw_p'
     *  Math: '<S307>/L^1//3'
     *  Product: '<S307>/Lug//V1'
     *  Product: '<S307>/w1'
     *  Product: '<S307>/w2'
     *  Sum: '<S307>/Sum'
     */
    pid_control_V1_B.w_o[0] = (pid_control_V1_B.chi / rt_powd_snf
      (pid_control_V1_B.Sum_hl, 0.33333333333333331) * pid_control_V1_ConstB.u16
      * pid_control_V1_B.Product[3] - pid_control_V1_X.pgw_p_CSTATE[0]) *
      pid_control_V1_B.sina;

    /* Math: '<S307>/L^1//3' */
    if (pid_control_V1_B.Sum5 < 0.0) {
      pid_control_V1_B.sinb = -rt_powd_snf(-pid_control_V1_B.Sum5,
        0.33333333333333331);
    } else {
      pid_control_V1_B.sinb = rt_powd_snf(pid_control_V1_B.Sum5,
        0.33333333333333331);
    }

    /* Product: '<S307>/w' incorporates:
     *  Fcn: '<S307>/sqrt(0.8//V)'
     *  Integrator: '<S307>/pgw_p'
     *  Math: '<S307>/L^1//3'
     *  Product: '<S307>/Lug//V1'
     *  Product: '<S307>/w1'
     *  Product: '<S307>/w2'
     *  Sum: '<S307>/Sum'
     */
    pid_control_V1_B.w_o[1] = (pid_control_V1_B.chi / pid_control_V1_B.sinb *
      pid_control_V1_ConstB.u16 * pid_control_V1_B.Product[3] -
      pid_control_V1_X.pgw_p_CSTATE[1]) * pid_control_V1_B.sina;

    /* Product: '<S307>/sigma_w' incorporates:
     *  Integrator: '<S307>/pgw_p'
     */
    tmp_2 = _mm_mul_pd(_mm_set_pd(pid_control_V1_B.Sum1_g,
      pid_control_V1_ConstB.sigma_wg), _mm_loadu_pd
                       (&pid_control_V1_X.pgw_p_CSTATE[0]));

    /* Product: '<S307>/sigma_w' */
    _mm_storeu_pd(&pid_control_V1_B.sigma_w[0], tmp_2);
  }

  /* End of Outputs for SubSystem: '<S296>/Hpgw' */

  /* Outputs for Enabled SubSystem: '<S297>/Hwgw(s)' incorporates:
   *  EnablePort: '<S312>/Enable'
   */
  if (tmp_0 && rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
      (!pid_control_V1_DW.Hwgws_MODE)) {
    (void) memset(&(pid_control_V1_XDis.wg_p1_CSTATE), 0,
                  4*sizeof(boolean_T));

    /* InitializeConditions for Integrator: '<S312>/wg_p1' */
    pid_control_V1_X.wg_p1_CSTATE[0] = 0.0;

    /* InitializeConditions for Integrator: '<S312>/wg_p2' */
    pid_control_V1_X.wg_p2_CSTATE[0] = 0.0;

    /* InitializeConditions for Integrator: '<S312>/wg_p1' */
    pid_control_V1_X.wg_p1_CSTATE[1] = 0.0;

    /* InitializeConditions for Integrator: '<S312>/wg_p2' */
    pid_control_V1_X.wg_p2_CSTATE[1] = 0.0;
    pid_control_V1_DW.Hwgws_MODE = true;
  }

  if (pid_control_V1_DW.Hwgws_MODE) {
    /* Product: '<S312>/Lwg//V' incorporates:
     *  Gain: '<S303>/Lw'
     */
    pid_control_V1_B.chi = pid_control_V1_B.Sum_hl /
      pid_control_V1_B.SignPreSat_a;

    /* Product: '<S312>/w' incorporates:
     *  Gain: '<S312>/1//pi'
     *  Integrator: '<S312>/wg_p1'
     *  Product: '<S312>/Lug//V1'
     *  Sqrt: '<S312>/sqrt1'
     *  Sum: '<S312>/Sum'
     */
    pid_control_V1_B.sina = (sqrt(0.31830988618379069 * pid_control_V1_B.chi) *
      pid_control_V1_B.Product[2] - pid_control_V1_X.wg_p1_CSTATE[0]) /
      pid_control_V1_B.chi;
    pid_control_V1_B.w[0] = pid_control_V1_B.sina;

    /* Product: '<S312>/w ' incorporates:
     *  Integrator: '<S312>/wg_p1'
     *  Integrator: '<S312>/wg_p2'
     *  Product: '<S312>/Lwg//V '
     *  Sum: '<S312>/Sum1'
     */
    pid_control_V1_B.w_a[0] = (pid_control_V1_B.sina *
      pid_control_V1_ConstB.sqrt_a * pid_control_V1_B.chi +
      (pid_control_V1_X.wg_p1_CSTATE[0] - pid_control_V1_X.wg_p2_CSTATE[0])) /
      pid_control_V1_B.chi;

    /* Product: '<S312>/Lwg//V' */
    pid_control_V1_B.chi = pid_control_V1_B.Sum5 / pid_control_V1_B.SignPreSat_a;

    /* Product: '<S312>/w' incorporates:
     *  Gain: '<S312>/1//pi'
     *  Integrator: '<S312>/wg_p1'
     *  Product: '<S312>/Lug//V1'
     *  Sqrt: '<S312>/sqrt1'
     *  Sum: '<S312>/Sum'
     */
    pid_control_V1_B.sina = (sqrt(0.31830988618379069 * pid_control_V1_B.chi) *
      pid_control_V1_B.Product[2] - pid_control_V1_X.wg_p1_CSTATE[1]) /
      pid_control_V1_B.chi;
    pid_control_V1_B.w[1] = pid_control_V1_B.sina;

    /* Product: '<S312>/w ' incorporates:
     *  Integrator: '<S312>/wg_p1'
     *  Integrator: '<S312>/wg_p2'
     *  Product: '<S312>/Lwg//V '
     *  Sum: '<S312>/Sum1'
     */
    pid_control_V1_B.w_a[1] = (pid_control_V1_B.sina *
      pid_control_V1_ConstB.sqrt_a * pid_control_V1_B.chi +
      (pid_control_V1_X.wg_p1_CSTATE[1] - pid_control_V1_X.wg_p2_CSTATE[1])) /
      pid_control_V1_B.chi;

    /* Product: '<S312>/Lwg//V 1' incorporates:
     *  Integrator: '<S312>/wg_p2'
     */
    tmp_2 = _mm_mul_pd(_mm_set_pd(pid_control_V1_B.Sum1_g,
      pid_control_V1_ConstB.sigma_wg), _mm_loadu_pd
                       (&pid_control_V1_X.wg_p2_CSTATE[0]));

    /* Product: '<S312>/Lwg//V 1' */
    _mm_storeu_pd(&pid_control_V1_B.LwgV1[0], tmp_2);
  }

  /* End of Outputs for SubSystem: '<S297>/Hwgw(s)' */

  /* Outputs for Enabled SubSystem: '<S296>/Hqgw' incorporates:
   *  EnablePort: '<S308>/Enable'
   */
  if (tmp_0 && rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
      (!pid_control_V1_DW.Hqgw_MODE)) {
    (void) memset(&(pid_control_V1_XDis.qgw_p_CSTATE), 0,
                  2*sizeof(boolean_T));

    /* InitializeConditions for Integrator: '<S308>/qgw_p' */
    pid_control_V1_X.qgw_p_CSTATE[0] = 0.0;
    pid_control_V1_X.qgw_p_CSTATE[1] = 0.0;
    pid_control_V1_DW.Hqgw_MODE = true;
  }

  if (pid_control_V1_DW.Hqgw_MODE) {
    /* Gain: '<S308>/pi//4' */
    pid_control_V1_B.chi = 0.78539816339744828 * pid_control_V1_B.SignPreSat_a;

    /* Product: '<S308>/w' incorporates:
     *  Integrator: '<S308>/qgw_p'
     *  Product: '<S308>/wg//V'
     *  Sum: '<S308>/Sum'
     */
    pid_control_V1_B.Sum5 = (pid_control_V1_B.LwgV1[0] /
      pid_control_V1_B.SignPreSat_a - pid_control_V1_X.qgw_p_CSTATE[0]) *
      (pid_control_V1_B.chi / pid_control_V1_ConstB.UnitConversion_n);
    pid_control_V1_B.w_e0[0] = pid_control_V1_B.Sum5;

    /* UnaryMinus: '<S308>/Unary Minus' */
    pid_control_V1_B.UnaryMinus[0] = -pid_control_V1_B.Sum5;

    /* Product: '<S308>/w' incorporates:
     *  Integrator: '<S308>/qgw_p'
     *  Product: '<S308>/wg//V'
     *  Sum: '<S308>/Sum'
     */
    pid_control_V1_B.Sum5 = (pid_control_V1_B.LwgV1[1] /
      pid_control_V1_B.SignPreSat_a - pid_control_V1_X.qgw_p_CSTATE[1]) *
      (pid_control_V1_B.chi / pid_control_V1_ConstB.UnitConversion_n);
    pid_control_V1_B.w_e0[1] = pid_control_V1_B.Sum5;

    /* UnaryMinus: '<S308>/Unary Minus' */
    pid_control_V1_B.UnaryMinus[1] = -pid_control_V1_B.Sum5;
  }

  /* End of Outputs for SubSystem: '<S296>/Hqgw' */

  /* Saturate: '<S314>/Limit Height h<1000ft' */
  if (pid_control_V1_B.SignPreSat > 1000.0) {
    pid_control_V1_B.chi = 1000.0;
  } else if (pid_control_V1_B.SignPreSat < 0.0) {
    pid_control_V1_B.chi = 0.0;
  } else {
    pid_control_V1_B.chi = pid_control_V1_B.SignPreSat;
  }

  /* Product: '<S314>/sigma_ug, sigma_vg' incorporates:
   *  Fcn: '<S314>/Low Altitude Intensity'
   *  Saturate: '<S314>/Limit Height h<1000ft'
   */
  pid_control_V1_B.sina = 1.0 / rt_powd_snf(0.000823 * pid_control_V1_B.chi +
    0.177, 0.4) * pid_control_V1_ConstB.sigma_wg;

  /* Fcn: '<S331>/Low Altitude Scale Length' */
  pid_control_V1_B.Sum_hl /= rt_powd_snf(0.000823 * pid_control_V1_B.Sum_hl +
    0.177, 1.2);

  /* Gain: '<S303>/Lv' */
  pid_control_V1_B.Sum5 = pid_control_V1_ConstB.UnitConversion_c;

  /* Outputs for Enabled SubSystem: '<S297>/Hvgw(s)' incorporates:
   *  EnablePort: '<S311>/Enable'
   */
  if (tmp_0 && rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
      (!pid_control_V1_DW.Hvgws_MODE)) {
    (void) memset(&(pid_control_V1_XDis.vg_p1_CSTATE), 0,
                  4*sizeof(boolean_T));

    /* InitializeConditions for Integrator: '<S311>/vg_p1' */
    pid_control_V1_X.vg_p1_CSTATE[0] = 0.0;

    /* InitializeConditions for Integrator: '<S311>/vgw_p2' */
    pid_control_V1_X.vgw_p2_CSTATE[0] = 0.0;

    /* InitializeConditions for Integrator: '<S311>/vg_p1' */
    pid_control_V1_X.vg_p1_CSTATE[1] = 0.0;

    /* InitializeConditions for Integrator: '<S311>/vgw_p2' */
    pid_control_V1_X.vgw_p2_CSTATE[1] = 0.0;
    pid_control_V1_DW.Hvgws_MODE = true;
  }

  if (pid_control_V1_DW.Hvgws_MODE) {
    /* Product: '<S311>/Lvg//V' incorporates:
     *  Gain: '<S303>/Lv'
     */
    pid_control_V1_B.chi = pid_control_V1_B.Sum_hl /
      pid_control_V1_B.SignPreSat_a;

    /* Product: '<S311>/w' incorporates:
     *  Gain: '<S311>/(1//pi)'
     *  Integrator: '<S311>/vg_p1'
     *  Product: '<S311>/Lug//V1'
     *  Sqrt: '<S311>/sqrt'
     *  Sum: '<S311>/Sum'
     */
    pid_control_V1_B.sinb = (sqrt(0.31830988618379069 * pid_control_V1_B.chi) *
      pid_control_V1_B.Product[1] - pid_control_V1_X.vg_p1_CSTATE[0]) /
      pid_control_V1_B.chi;
    pid_control_V1_B.w_g[0] = pid_control_V1_B.sinb;

    /* Product: '<S311>/w ' incorporates:
     *  Gain: '<S311>/sqrt(3)'
     *  Integrator: '<S311>/vg_p1'
     *  Integrator: '<S311>/vgw_p2'
     *  Product: '<S311>/Lvg//V '
     *  Sum: '<S311>/Sum1'
     */
    pid_control_V1_B.w_e[0] = (pid_control_V1_B.sinb * pid_control_V1_B.chi *
      1.7320508075688772 + (pid_control_V1_X.vg_p1_CSTATE[0] -
      pid_control_V1_X.vgw_p2_CSTATE[0])) / pid_control_V1_B.chi;

    /* Product: '<S311>/Lvg//V' */
    pid_control_V1_B.chi = pid_control_V1_B.Sum5 / pid_control_V1_B.SignPreSat_a;

    /* Product: '<S311>/w' incorporates:
     *  Gain: '<S311>/(1//pi)'
     *  Integrator: '<S311>/vg_p1'
     *  Product: '<S311>/Lug//V1'
     *  Sqrt: '<S311>/sqrt'
     *  Sum: '<S311>/Sum'
     */
    pid_control_V1_B.sinb = (sqrt(0.31830988618379069 * pid_control_V1_B.chi) *
      pid_control_V1_B.Product[1] - pid_control_V1_X.vg_p1_CSTATE[1]) /
      pid_control_V1_B.chi;
    pid_control_V1_B.w_g[1] = pid_control_V1_B.sinb;

    /* Product: '<S311>/w ' incorporates:
     *  Gain: '<S311>/sqrt(3)'
     *  Integrator: '<S311>/vg_p1'
     *  Integrator: '<S311>/vgw_p2'
     *  Product: '<S311>/Lvg//V '
     *  Sum: '<S311>/Sum1'
     */
    pid_control_V1_B.w_e[1] = (pid_control_V1_B.sinb * pid_control_V1_B.chi *
      1.7320508075688772 + (pid_control_V1_X.vg_p1_CSTATE[1] -
      pid_control_V1_X.vgw_p2_CSTATE[1])) / pid_control_V1_B.chi;

    /* Product: '<S311>/w 1' incorporates:
     *  Integrator: '<S311>/vgw_p2'
     */
    tmp_2 = _mm_mul_pd(_mm_set_pd(pid_control_V1_B.Sum1_g, pid_control_V1_B.sina),
                       _mm_loadu_pd(&pid_control_V1_X.vgw_p2_CSTATE[0]));

    /* Product: '<S311>/w 1' */
    _mm_storeu_pd(&pid_control_V1_B.w1[0], tmp_2);
  }

  /* End of Outputs for SubSystem: '<S297>/Hvgw(s)' */

  /* Outputs for Enabled SubSystem: '<S296>/Hrgw' incorporates:
   *  EnablePort: '<S309>/Enable'
   */
  if (tmp_0 && rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
      (!pid_control_V1_DW.Hrgw_MODE)) {
    (void) memset(&(pid_control_V1_XDis.rgw_p_CSTATE), 0,
                  2*sizeof(boolean_T));

    /* InitializeConditions for Integrator: '<S309>/rgw_p' */
    pid_control_V1_X.rgw_p_CSTATE[0] = 0.0;
    pid_control_V1_X.rgw_p_CSTATE[1] = 0.0;
    pid_control_V1_DW.Hrgw_MODE = true;
  }

  if (pid_control_V1_DW.Hrgw_MODE) {
    /* Product: '<S309>/vg//V' incorporates:
     *  Gain: '<S309>/pi//3'
     *  Integrator: '<S309>/rgw_p'
     *  Product: '<S309>/w'
     */
    tmp_2 = _mm_mul_pd(_mm_sub_pd(_mm_div_pd(_mm_loadu_pd(&pid_control_V1_B.w1[0]),
      _mm_set1_pd(pid_control_V1_B.SignPreSat_a)), _mm_loadu_pd
      (&pid_control_V1_X.rgw_p_CSTATE[0])), _mm_div_pd(_mm_set1_pd
      (1.0471975511965976 * pid_control_V1_B.SignPreSat_a), _mm_set1_pd
      (pid_control_V1_ConstB.UnitConversion_n)));

    /* Product: '<S309>/w' */
    _mm_storeu_pd(&pid_control_V1_B.w_d[0], tmp_2);
  }

  /* End of Outputs for SubSystem: '<S296>/Hrgw' */

  /* Outputs for Enabled SubSystem: '<S297>/Hugw(s)' incorporates:
   *  EnablePort: '<S310>/Enable'
   */
  if (tmp_0 && rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)->solverInfo) &&
      (!pid_control_V1_DW.Hugws_MODE)) {
    (void) memset(&(pid_control_V1_XDis.ug_p_CSTATE), 0,
                  2*sizeof(boolean_T));

    /* InitializeConditions for Integrator: '<S310>/ug_p' */
    pid_control_V1_X.ug_p_CSTATE[0] = 0.0;
    pid_control_V1_X.ug_p_CSTATE[1] = 0.0;
    pid_control_V1_DW.Hugws_MODE = true;
  }

  if (pid_control_V1_DW.Hugws_MODE) {
    /* Product: '<S310>/Lug//V' */
    pid_control_V1_B.chi = pid_control_V1_B.Sum_hl /
      pid_control_V1_B.SignPreSat_a;
    pid_control_V1_B.Sum5 = pid_control_V1_ConstB.UnitConversion_c /
      pid_control_V1_B.SignPreSat_a;

    /* Sqrt: '<S310>/sqrt' incorporates:
     *  Gain: '<S310>/(2//pi)'
     *  Integrator: '<S310>/ug_p'
     *  Product: '<S310>/Lug//V1'
     */
    tmp_2 = _mm_div_pd(_mm_sub_pd(_mm_mul_pd(_mm_set_pd(sqrt(0.63661977236758138
      * pid_control_V1_B.Sum5), sqrt(0.63661977236758138 * pid_control_V1_B.chi)),
      _mm_set1_pd(pid_control_V1_B.Product[0])), _mm_loadu_pd
      (&pid_control_V1_X.ug_p_CSTATE[0])), _mm_set_pd(pid_control_V1_B.Sum5,
      pid_control_V1_B.chi));

    /* Product: '<S310>/w' */
    _mm_storeu_pd(&pid_control_V1_B.w_n[0], tmp_2);

    /* Integrator: '<S310>/ug_p' incorporates:
     *  Product: '<S310>/w1'
     */
    tmp_2 = _mm_mul_pd(_mm_loadu_pd(&pid_control_V1_X.ug_p_CSTATE[0]),
                       _mm_set_pd(pid_control_V1_B.Sum1_g, pid_control_V1_B.sina));

    /* Product: '<S310>/w1' */
    _mm_storeu_pd(&pid_control_V1_B.w1_c[0], tmp_2);
  }

  /* End of Outputs for SubSystem: '<S297>/Hugw(s)' */

  /* Angle2Dcm: '<S16>/Rotation Angles to Direction Cosine Matrix' */
  pid_control_V1_B.chi = cos(pid_control_V1_B.x[6]);
  pid_control_V1_B.Sum1_g = sin(pid_control_V1_B.x[6]);
  pid_control_V1_B.Sum_hl = -sin(pid_control_V1_B.x[6]);
  pid_control_V1_B.Sum5 = cos(pid_control_V1_B.x[6]);
  pid_control_V1_B.Dtot = cos(pid_control_V1_B.x[7]);
  pid_control_V1_B.u2 = -sin(pid_control_V1_B.x[7]);
  pid_control_V1_B.SignPreSat_a = sin(pid_control_V1_B.x[7]);
  pid_control_V1_B.sina = cos(pid_control_V1_B.x[7]);
  pid_control_V1_B.sinb = cos(pid_control_V1_B.x[8]);
  pid_control_V1_B.sinc = sin(pid_control_V1_B.x[8]);
  pid_control_V1_B.cosa = -sin(pid_control_V1_B.x[8]);
  pid_control_V1_B.Ltot = cos(pid_control_V1_B.x[8]);
  pid_control_V1_B.cosb = 0.0 * pid_control_V1_B.SignPreSat_a +
    pid_control_V1_B.Dtot;
  pid_control_V1_B.cosc = 0.0 * pid_control_V1_B.sina + pid_control_V1_B.u2;
  pid_control_V1_B.CQ = pid_control_V1_B.sinb * 0.0;
  pid_control_V1_B.Cl = 0.0 * pid_control_V1_B.Dtot;
  pid_control_V1_B.Dtot = (pid_control_V1_B.Cl + pid_control_V1_B.CQ) +
    pid_control_V1_B.sinc * pid_control_V1_B.SignPreSat_a;
  pid_control_V1_B.sinb += pid_control_V1_B.sinc * 0.0;
  pid_control_V1_B.u2 *= 0.0;
  pid_control_V1_B.sinc = (pid_control_V1_B.u2 + pid_control_V1_B.CQ) +
    pid_control_V1_B.sinc * pid_control_V1_B.sina;
  pid_control_V1_B.CQ = pid_control_V1_B.cosa * 0.0;
  pid_control_V1_B.SignPreSat_a = (pid_control_V1_B.Cl + pid_control_V1_B.CQ) +
    pid_control_V1_B.SignPreSat_a * pid_control_V1_B.Ltot;
  pid_control_V1_B.cosa += pid_control_V1_B.Ltot * 0.0;
  pid_control_V1_B.sina = (pid_control_V1_B.u2 + pid_control_V1_B.CQ) +
    pid_control_V1_B.Ltot * pid_control_V1_B.sina;
  pid_control_V1_B.Ltot = pid_control_V1_B.cosc * 0.0;
  pid_control_V1_B.RotationAnglestoDirectionCo[0] = (pid_control_V1_B.cosb *
    pid_control_V1_B.chi + 0.0 * pid_control_V1_B.Sum_hl) +
    pid_control_V1_B.Ltot;
  pid_control_V1_B.CQ = pid_control_V1_B.sinc * 0.0;
  pid_control_V1_B.RotationAnglestoDirectionCo[1] = (pid_control_V1_B.chi *
    pid_control_V1_B.Dtot + pid_control_V1_B.Sum_hl * pid_control_V1_B.sinb) +
    pid_control_V1_B.CQ;
  pid_control_V1_B.Cl = pid_control_V1_B.sina * 0.0;
  pid_control_V1_B.RotationAnglestoDirectionCo[2] = (pid_control_V1_B.chi *
    pid_control_V1_B.SignPreSat_a + pid_control_V1_B.Sum_hl *
    pid_control_V1_B.cosa) + pid_control_V1_B.Cl;
  pid_control_V1_B.RotationAnglestoDirectionCo[3] = (pid_control_V1_B.cosb *
    pid_control_V1_B.Sum1_g + 0.0 * pid_control_V1_B.Sum5) +
    pid_control_V1_B.Ltot;
  pid_control_V1_B.RotationAnglestoDirectionCo[4] = (pid_control_V1_B.Sum1_g *
    pid_control_V1_B.Dtot + pid_control_V1_B.sinb * pid_control_V1_B.Sum5) +
    pid_control_V1_B.CQ;
  pid_control_V1_B.RotationAnglestoDirectionCo[5] = (pid_control_V1_B.Sum1_g *
    pid_control_V1_B.SignPreSat_a + pid_control_V1_B.Sum5 *
    pid_control_V1_B.cosa) + pid_control_V1_B.Cl;
  pid_control_V1_B.RotationAnglestoDirectionCo[6] = pid_control_V1_B.cosb * 0.0
    + pid_control_V1_B.cosc;
  pid_control_V1_B.RotationAnglestoDirectionCo[7] = (pid_control_V1_B.Dtot * 0.0
    + pid_control_V1_B.sinb * 0.0) + pid_control_V1_B.sinc;
  pid_control_V1_B.RotationAnglestoDirectionCo[8] =
    (pid_control_V1_B.SignPreSat_a * 0.0 + pid_control_V1_B.cosa * 0.0) +
    pid_control_V1_B.sina;

  /* If: '<S301>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' incorporates:
   *  Constant: '<S315>/max_height_low'
   *  If: '<S302>/if Height < Max low altitude  elseif Height > Min isotropic altitude '
   *  Product: '<S315>/Product1'
   *  Product: '<S320>/Product1'
   *  Product: '<S320>/Product2'
   *  Product: '<S322>/Product1'
   *  Product: '<S322>/Product2'
   *  Sum: '<S315>/Sum1'
   *  Sum: '<S315>/Sum2'
   *  Sum: '<S315>/Sum3'
   *  Sum: '<S320>/Sum'
   *  Sum: '<S322>/Sum'
   */
  rtPrevAction = pid_control_V1_DW.ifHeightMaxlowaltitudeelseifHei;
  serverAvailableOnTime = rtsiIsModeUpdateTimeStep(&(&pid_control_V1_M)
    ->solverInfo);
  if (serverAvailableOnTime) {
    if (pid_control_V1_B.SignPreSat <= 1000.0) {
      rtAction = 0;
    } else if (pid_control_V1_B.SignPreSat >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    pid_control_V1_DW.ifHeightMaxlowaltitudeelseifHei = rtAction;
  } else {
    rtAction = pid_control_V1_DW.ifHeightMaxlowaltitudeelseifHei;
  }

  if (rtPrevAction != rtAction) {
    rtsiSetBlockStateForSolverChangedAtMajorStep(&(&pid_control_V1_M)
      ->solverInfo, true);
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S301>/Low altitude  rates' incorporates:
     *  ActionPort: '<S316>/Action Port'
     */
    /* SignalConversion generated from: '<S321>/Vector Concatenate' */
    pid_control_V1_B.Product_be[2] = pid_control_V1_B.w_d[0];

    /* Trigonometry: '<S322>/Trigonometric Function1' incorporates:
     *  UnitConversion: '<S295>/Unit Conversion'
     */
    pid_control_V1_B.chi = sin(pid_control_V1_ConstB.UnitConversion);
    pid_control_V1_B.Sum1_g = cos(pid_control_V1_ConstB.UnitConversion);
    _mm_storeu_pd(&pid_control_V1_B.Product_be[0], _mm_add_pd(_mm_mul_pd
      (_mm_set_pd(pid_control_V1_B.chi, pid_control_V1_B.sigma_w[0]), _mm_set_pd
       (pid_control_V1_B.sigma_w[0], pid_control_V1_B.Sum1_g)), _mm_mul_pd
      (_mm_mul_pd(_mm_set_pd(pid_control_V1_B.UnaryMinus[0],
      pid_control_V1_B.chi), _mm_set_pd(pid_control_V1_B.Sum1_g,
      pid_control_V1_B.UnaryMinus[0])), _mm_set_pd(1.0, -1.0))));

    /* Product: '<S321>/Product' incorporates:
     *  Angle2Dcm: '<S16>/Rotation Angles to Direction Cosine Matrix'
     *  Concatenate: '<S321>/Vector Concatenate'
     *  Product: '<S322>/Product1'
     *  Product: '<S322>/Product2'
     *  Reshape: '<S321>/Reshape1'
     *  Sum: '<S322>/Sum'
     */
    pid_control_V1_B.chi = 0.0;
    pid_control_V1_B.Sum1_g = 0.0;
    pid_control_V1_B.Sum_hl = 0.0;
    for (i = 0; i < 3; i++) {
      tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&pid_control_V1_B.RotationAnglestoDirectionCo[3 * i]), _mm_set1_pd
        (pid_control_V1_B.Product_be[i])), _mm_set_pd(pid_control_V1_B.Sum1_g,
        pid_control_V1_B.chi));
      _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
      pid_control_V1_B.chi = pid_control_V1_B.dv[0];
      pid_control_V1_B.Sum1_g = pid_control_V1_B.dv[1];
      pid_control_V1_B.Sum_hl += pid_control_V1_B.RotationAnglestoDirectionCo[3 *
        i + 2] * pid_control_V1_B.Product_be[i];
    }

    pid_control_V1_B.wbe_b[2] = pid_control_V1_B.Sum_hl;
    pid_control_V1_B.wbe_b[1] = pid_control_V1_B.Sum1_g;
    pid_control_V1_B.wbe_b[0] = pid_control_V1_B.chi;

    /* End of Product: '<S321>/Product' */
    /* End of Outputs for SubSystem: '<S301>/Low altitude  rates' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S301>/Medium//High  altitude rates' incorporates:
     *  ActionPort: '<S317>/Action Port'
     */
    /* Gain: '<S317>/Gain' */
    pid_control_V1_B.wbe_b[0] = pid_control_V1_B.sigma_w[1];
    pid_control_V1_B.wbe_b[1] = pid_control_V1_B.UnaryMinus[1];
    pid_control_V1_B.wbe_b[2] = pid_control_V1_B.w_d[1];

    /* End of Outputs for SubSystem: '<S301>/Medium//High  altitude rates' */
    break;

   default:
    /* Outputs for IfAction SubSystem: '<S301>/Interpolate  rates' incorporates:
     *  ActionPort: '<S315>/Action Port'
     */
    /* Trigonometry: '<S320>/Trigonometric Function' incorporates:
     *  UnitConversion: '<S295>/Unit Conversion'
     */
    pid_control_V1_B.chi = sin(pid_control_V1_ConstB.UnitConversion);
    pid_control_V1_B.Sum1_g = cos(pid_control_V1_ConstB.UnitConversion);
    _mm_storeu_pd(&pid_control_V1_B.wbe_b[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (pid_control_V1_B.chi, pid_control_V1_B.sigma_w[0]), _mm_set_pd
      (pid_control_V1_B.sigma_w[0], pid_control_V1_B.Sum1_g)), _mm_mul_pd
      (_mm_mul_pd(_mm_set_pd(pid_control_V1_B.UnaryMinus[0],
      pid_control_V1_B.chi), _mm_set_pd(pid_control_V1_B.Sum1_g,
      pid_control_V1_B.UnaryMinus[0])), _mm_set_pd(1.0, -1.0))));

    /* SignalConversion generated from: '<S319>/Vector Concatenate' incorporates:
     *  Product: '<S320>/Product1'
     *  Product: '<S320>/Product2'
     *  Sum: '<S320>/Sum'
     */
    pid_control_V1_B.wbe_b[2] = pid_control_V1_B.w_d[0];

    /* Product: '<S319>/Product' incorporates:
     *  Angle2Dcm: '<S16>/Rotation Angles to Direction Cosine Matrix'
     *  Concatenate: '<S319>/Vector Concatenate'
     */
    pid_control_V1_B.chi = 0.0;
    pid_control_V1_B.Sum1_g = 0.0;
    pid_control_V1_B.Sum_hl = 0.0;
    for (i = 0; i < 3; i++) {
      tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&pid_control_V1_B.RotationAnglestoDirectionCo[3 * i]), _mm_set1_pd
        (pid_control_V1_B.wbe_b[i])), _mm_set_pd(pid_control_V1_B.Sum1_g,
        pid_control_V1_B.chi));
      _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
      pid_control_V1_B.chi = pid_control_V1_B.dv[0];
      pid_control_V1_B.Sum1_g = pid_control_V1_B.dv[1];
      pid_control_V1_B.Sum_hl += pid_control_V1_B.RotationAnglestoDirectionCo[3 *
        i + 2] * pid_control_V1_B.wbe_b[i];
    }

    pid_control_V1_B.Product_be[2] = pid_control_V1_B.Sum_hl;
    pid_control_V1_B.Product_be[1] = pid_control_V1_B.Sum1_g;
    pid_control_V1_B.Product_be[0] = pid_control_V1_B.chi;
    tmp_2 = _mm_add_pd(_mm_div_pd(_mm_mul_pd(_mm_sub_pd(_mm_set_pd
      (pid_control_V1_B.UnaryMinus[1], pid_control_V1_B.sigma_w[1]),
      _mm_loadu_pd(&pid_control_V1_B.Product_be[0])), _mm_sub_pd(_mm_set1_pd
      (pid_control_V1_B.SignPreSat), _mm_set1_pd(1000.0))), _mm_set1_pd
      (pid_control_V1_ConstB.Sum_a)), _mm_loadu_pd(&pid_control_V1_B.Product_be
      [0]));
    _mm_storeu_pd(&pid_control_V1_B.wbe_b[0], tmp_2);

    /* Sum: '<S315>/Sum3' incorporates:
     *  Constant: '<S315>/max_height_low'
     *  Product: '<S315>/Product1'
     *  Product: '<S319>/Product'
     *  Sum: '<S315>/Sum1'
     *  Sum: '<S315>/Sum2'
     */
    pid_control_V1_B.wbe_b[2] = (pid_control_V1_B.w_d[1] -
      pid_control_V1_B.Sum_hl) * (pid_control_V1_B.SignPreSat - 1000.0) /
      pid_control_V1_ConstB.Sum_a + pid_control_V1_B.Sum_hl;

    /* End of Outputs for SubSystem: '<S301>/Interpolate  rates' */
    break;
  }

  /* End of If: '<S301>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */

  /* If: '<S302>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' incorporates:
   *  Constant: '<S323>/max_height_low'
   *  Product: '<S323>/Product1'
   *  Product: '<S328>/Product1'
   *  Product: '<S328>/Product2'
   *  Product: '<S330>/Product1'
   *  Product: '<S330>/Product2'
   *  Sum: '<S323>/Sum1'
   *  Sum: '<S323>/Sum2'
   *  Sum: '<S323>/Sum3'
   *  Sum: '<S328>/Sum'
   *  Sum: '<S330>/Sum'
   */
  rtPrevAction = pid_control_V1_DW.ifHeightMaxlowaltitudeelseifH_k;
  if (serverAvailableOnTime) {
    if (pid_control_V1_B.SignPreSat <= 1000.0) {
      rtAction = 0;
    } else if (pid_control_V1_B.SignPreSat >= 2000.0) {
      rtAction = 1;
    } else {
      rtAction = 2;
    }

    pid_control_V1_DW.ifHeightMaxlowaltitudeelseifH_k = rtAction;
  } else {
    rtAction = pid_control_V1_DW.ifHeightMaxlowaltitudeelseifH_k;
  }

  if (rtPrevAction != rtAction) {
    rtsiSetBlockStateForSolverChangedAtMajorStep(&(&pid_control_V1_M)
      ->solverInfo, true);
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S302>/Low altitude  velocities' incorporates:
     *  ActionPort: '<S324>/Action Port'
     */
    /* SignalConversion generated from: '<S329>/Vector Concatenate' */
    pid_control_V1_B.FA_b[2] = pid_control_V1_B.LwgV1[0];

    /* Trigonometry: '<S330>/Trigonometric Function' incorporates:
     *  UnitConversion: '<S295>/Unit Conversion'
     */
    pid_control_V1_B.chi = sin(pid_control_V1_ConstB.UnitConversion);
    pid_control_V1_B.SignPreSat = cos(pid_control_V1_ConstB.UnitConversion);
    _mm_storeu_pd(&pid_control_V1_B.FA_b[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
      (pid_control_V1_B.chi, pid_control_V1_B.w1_c[0]), _mm_set_pd
      (pid_control_V1_B.w1_c[0], pid_control_V1_B.SignPreSat)), _mm_mul_pd
      (_mm_mul_pd(_mm_set_pd(pid_control_V1_B.w1[0], pid_control_V1_B.chi),
                  _mm_set_pd(pid_control_V1_B.SignPreSat, pid_control_V1_B.w1[0])),
       _mm_set_pd(1.0, -1.0))));

    /* Product: '<S329>/Product' incorporates:
     *  Angle2Dcm: '<S16>/Rotation Angles to Direction Cosine Matrix'
     *  Concatenate: '<S329>/Vector Concatenate'
     *  Product: '<S330>/Product1'
     *  Product: '<S330>/Product2'
     *  Reshape: '<S329>/Reshape1'
     *  Sum: '<S330>/Sum'
     */
    pid_control_V1_B.chi = 0.0;
    pid_control_V1_B.Sum1_g = 0.0;
    pid_control_V1_B.Sum_hl = 0.0;
    for (i = 0; i < 3; i++) {
      tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&pid_control_V1_B.RotationAnglestoDirectionCo[3 * i]), _mm_set1_pd
        (pid_control_V1_B.FA_b[i])), _mm_set_pd(pid_control_V1_B.Sum1_g,
        pid_control_V1_B.chi));
      _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
      pid_control_V1_B.chi = pid_control_V1_B.dv[0];
      pid_control_V1_B.Sum1_g = pid_control_V1_B.dv[1];
      pid_control_V1_B.Sum_hl += pid_control_V1_B.RotationAnglestoDirectionCo[3 *
        i + 2] * pid_control_V1_B.FA_b[i];
    }

    pid_control_V1_B.Product_be[2] = pid_control_V1_B.Sum_hl;
    pid_control_V1_B.Product_be[1] = pid_control_V1_B.Sum1_g;
    pid_control_V1_B.Product_be[0] = pid_control_V1_B.chi;

    /* End of Product: '<S329>/Product' */
    /* End of Outputs for SubSystem: '<S302>/Low altitude  velocities' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S302>/Medium//High  altitude velocities' incorporates:
     *  ActionPort: '<S325>/Action Port'
     */
    /* Gain: '<S325>/Gain' */
    pid_control_V1_B.Product_be[0] = pid_control_V1_B.w1_c[1];
    pid_control_V1_B.Product_be[1] = pid_control_V1_B.w1[1];
    pid_control_V1_B.Product_be[2] = pid_control_V1_B.LwgV1[1];

    /* End of Outputs for SubSystem: '<S302>/Medium//High  altitude velocities' */
    break;

   default:
    /* Outputs for IfAction SubSystem: '<S302>/Interpolate  velocities' incorporates:
     *  ActionPort: '<S323>/Action Port'
     */
    /* Trigonometry: '<S328>/Trigonometric Function' incorporates:
     *  UnitConversion: '<S295>/Unit Conversion'
     */
    pid_control_V1_B.chi = sin(pid_control_V1_ConstB.UnitConversion);
    pid_control_V1_B.Sum1_g = cos(pid_control_V1_ConstB.UnitConversion);
    _mm_storeu_pd(&pid_control_V1_B.Product_be[0], _mm_add_pd(_mm_mul_pd
      (_mm_set_pd(pid_control_V1_B.chi, pid_control_V1_B.w1_c[0]), _mm_set_pd
       (pid_control_V1_B.w1_c[0], pid_control_V1_B.Sum1_g)), _mm_mul_pd
      (_mm_mul_pd(_mm_set_pd(pid_control_V1_B.w1[0], pid_control_V1_B.chi),
                  _mm_set_pd(pid_control_V1_B.Sum1_g, pid_control_V1_B.w1[0])),
       _mm_set_pd(1.0, -1.0))));

    /* SignalConversion generated from: '<S327>/Vector Concatenate' incorporates:
     *  Product: '<S328>/Product1'
     *  Product: '<S328>/Product2'
     *  Sum: '<S328>/Sum'
     */
    pid_control_V1_B.Product_be[2] = pid_control_V1_B.LwgV1[0];

    /* Product: '<S327>/Product' incorporates:
     *  Angle2Dcm: '<S16>/Rotation Angles to Direction Cosine Matrix'
     *  Concatenate: '<S327>/Vector Concatenate'
     */
    pid_control_V1_B.Sum1_g = 0.0;
    pid_control_V1_B.Sum_hl = 0.0;
    pid_control_V1_B.Sum5 = 0.0;
    for (i = 0; i < 3; i++) {
      tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&pid_control_V1_B.RotationAnglestoDirectionCo[3 * i]), _mm_set1_pd
        (pid_control_V1_B.Product_be[i])), _mm_set_pd(pid_control_V1_B.Sum_hl,
        pid_control_V1_B.Sum1_g));
      _mm_storeu_pd(&pid_control_V1_B.dv[0], tmp_2);
      pid_control_V1_B.Sum1_g = pid_control_V1_B.dv[0];
      pid_control_V1_B.Sum_hl = pid_control_V1_B.dv[1];
      pid_control_V1_B.Sum5 += pid_control_V1_B.RotationAnglestoDirectionCo[3 *
        i + 2] * pid_control_V1_B.Product_be[i];
    }

    pid_control_V1_B.FA_b[2] = pid_control_V1_B.Sum5;
    pid_control_V1_B.FA_b[1] = pid_control_V1_B.Sum_hl;
    pid_control_V1_B.FA_b[0] = pid_control_V1_B.Sum1_g;
    tmp_2 = _mm_add_pd(_mm_div_pd(_mm_mul_pd(_mm_sub_pd(_mm_set_pd
      (pid_control_V1_B.w1[1], pid_control_V1_B.w1_c[1]), _mm_loadu_pd
      (&pid_control_V1_B.FA_b[0])), _mm_sub_pd(_mm_set1_pd
      (pid_control_V1_B.SignPreSat), _mm_set1_pd(1000.0))), _mm_set1_pd
      (pid_control_V1_ConstB.Sum)), _mm_loadu_pd(&pid_control_V1_B.FA_b[0]));
    _mm_storeu_pd(&pid_control_V1_B.Product_be[0], tmp_2);

    /* Sum: '<S323>/Sum3' incorporates:
     *  Constant: '<S323>/max_height_low'
     *  Product: '<S323>/Product1'
     *  Product: '<S327>/Product'
     *  Sum: '<S323>/Sum1'
     *  Sum: '<S323>/Sum2'
     */
    pid_control_V1_B.Product_be[2] = (pid_control_V1_B.LwgV1[1] -
      pid_control_V1_B.Sum5) * (pid_control_V1_B.SignPreSat - 1000.0) /
      pid_control_V1_ConstB.Sum + pid_control_V1_B.Sum5;

    /* End of Outputs for SubSystem: '<S302>/Interpolate  velocities' */
    break;
  }

  /* UnitConversion: '<S289>/Unit Conversion' */
  /* Unit Conversion - from: ft/s to: m/s
     Expression: output = (0.3048*input) + (0) */
  tmp_2 = _mm_mul_pd(_mm_set1_pd(0.3048), _mm_loadu_pd
                     (&pid_control_V1_B.Product_be[0]));
  _mm_storeu_pd(&pid_control_V1_B.Product_be[0], tmp_2);
  pid_control_V1_B.Product_be[2] *= 0.3048;
  if (tmp_0) {
    /* MATLABSystem: '<S291>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_h = Sub_pid_control_V1_417.getLatestMessage(
      &rtb_SourceBlock_o2_j);

    /* Outputs for Enabled SubSystem: '<S291>/Enabled Subsystem' */
    pid_control__EnabledSubsystem_p(pid_control_V1_B.SourceBlock_o1_h,
      &rtb_SourceBlock_o2_j, &pid_control_V1_B.EnabledSubsystem_p);

    /* End of Outputs for SubSystem: '<S291>/Enabled Subsystem' */

    /* MATLABSystem: '<S292>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_k = Sub_pid_control_V1_423.getLatestMessage(
      &rtb_SourceBlock_o2_dd);

    /* Outputs for Enabled SubSystem: '<S292>/Enabled Subsystem' */
    pid_control__EnabledSubsystem_p(pid_control_V1_B.SourceBlock_o1_k,
      &rtb_SourceBlock_o2_dd, &pid_control_V1_B.EnabledSubsystem_g);

    /* End of Outputs for SubSystem: '<S292>/Enabled Subsystem' */

    /* MATLABSystem: '<S293>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1_c = Sub_pid_control_V1_443.getLatestMessage(
      &pid_control_V1_B.SourceBlock_o2_p);

    /* Outputs for Enabled SubSystem: '<S293>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1_c,
      &pid_control_V1_B.SourceBlock_o2_p, &pid_control_V1_B.EnabledSubsystem_k);

    /* End of Outputs for SubSystem: '<S293>/Enabled Subsystem' */

    /* MATLABSystem: '<S294>/SourceBlock' */
    pid_control_V1_B.SourceBlock_o1 = Sub_pid_control_V1_445.getLatestMessage
      (&pid_control_V1_B.SourceBlock_o2_k);

    /* Outputs for Enabled SubSystem: '<S294>/Enabled Subsystem' */
    pid_control_V1_EnabledSubsystem(pid_control_V1_B.SourceBlock_o1,
      &pid_control_V1_B.SourceBlock_o2_k, &pid_control_V1_B.EnabledSubsystem_pu);

    /* End of Outputs for SubSystem: '<S294>/Enabled Subsystem' */
  }

  /* Switch: '<S16>/Switch' incorporates:
   *  Constant: '<S16>/Constant'
   */
  if (!pid_control_V1_B.EnabledSubsystem_p.In1.data) {
    pid_control_V1_B.Product_be[0] = 0.0;
    pid_control_V1_B.Product_be[1] = 0.0;
    pid_control_V1_B.Product_be[2] = 0.0;
  }

  /* End of Switch: '<S16>/Switch' */

  /* TransferFcn: '<S16>/Transfer Fcn' */
  pid_control_V1_B.Sum_hl = 0.5303 * pid_control_V1_X.TransferFcn_CSTATE[0] +
    0.0 * pid_control_V1_X.TransferFcn_CSTATE[1];

  /* Switch: '<S16>/Switch2' incorporates:
   *  Constant: '<S16>/Constant3'
   */
  if (!(pid_control_V1_B.EnabledSubsystem_k.In1.data != 0.0)) {
    pid_control_V1_B.Sum_hl = 0.0;
  }

  /* End of Switch: '<S16>/Switch2' */

  /* Sum: '<S16>/Sum' */
  pid_control_V1_B.Sum[0] = pid_control_V1_B.Product_be[0];
  pid_control_V1_B.Sum[1] = pid_control_V1_B.Product_be[1] +
    pid_control_V1_B.Sum_hl;
  pid_control_V1_B.Sum[2] = pid_control_V1_B.Product_be[2];

  /* Switch: '<S16>/Switch1' incorporates:
   *  Constant: '<S16>/Constant2'
   */
  if (!pid_control_V1_B.EnabledSubsystem_g.In1.data) {
    pid_control_V1_B.wbe_b[0] = 0.0;
    pid_control_V1_B.wbe_b[1] = 0.0;
    pid_control_V1_B.wbe_b[2] = 0.0;
  }

  /* End of Switch: '<S16>/Switch1' */

  /* Sum: '<S16>/Sum1' */
  pid_control_V1_B.Sum1[0] = pid_control_V1_B.wbe_b[0];

  /* Switch: '<S16>/Switch3' incorporates:
   *  Constant: '<S16>/Constant4'
   *  TransferFcn: '<S16>/Transfer Fcn1'
   */
  if (pid_control_V1_B.EnabledSubsystem_pu.In1.data != 0.0) {
    pid_control_V1_B.chi = -0.0003571 * pid_control_V1_X.TransferFcn1_CSTATE +
      0.03571 * pid_control_V1_B.Output;
  } else {
    pid_control_V1_B.chi = 0.0;
  }

  /* Sum: '<S16>/Sum1' incorporates:
   *  Switch: '<S16>/Switch3'
   */
  pid_control_V1_B.Sum1[1] = pid_control_V1_B.wbe_b[1] + pid_control_V1_B.chi;
  pid_control_V1_B.Sum1[2] = pid_control_V1_B.wbe_b[2];
  if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
    if (rtmIsMajorTimeStep((&pid_control_V1_M))) {
      /* Update for UnitDelay: '<Root>/Unit Delay3' */
      pid_control_V1_DW.UnitDelay3_DSTATE = pid_control_V1_B.Switch3;

      /* Update for UnitDelay: '<Root>/Unit Delay4' */
      pid_control_V1_DW.UnitDelay4_DSTATE = pid_control_V1_B.Switch4;

      /* Update for UnitDelay: '<Root>/Unit Delay5' */
      pid_control_V1_DW.UnitDelay5_DSTATE = pid_control_V1_B.Switch5;

      /* Update for UnitDelay: '<Root>/Unit Delay6' */
      pid_control_V1_DW.UnitDelay6_DSTATE = pid_control_V1_B.Switch6;

      /* Update for UnitDelay: '<Root>/Unit Delay2' */
      pid_control_V1_DW.UnitDelay2_DSTATE = pid_control_V1_B.Switch2;

      /* Update for Memory: '<S43>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput_o = pid_control_V1_B.AND3;

      /* Update for Memory: '<S97>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput_m = pid_control_V1_B.AND3_e;

      /* Update for Memory: '<S255>/Memory' */
      pid_control_V1_DW.Memory_PreviousInput_a = pid_control_V1_B.AND3_c;

      /* Update for Memory: '<S16>/Memory' incorporates:
       *  Sum: '<S16>/Sum'
       */
      pid_control_V1_DW.Memory_PreviousInput[0] = pid_control_V1_B.Sum[0];

      /* Update for Memory: '<S16>/Memory1' incorporates:
       *  Sum: '<S16>/Sum1'
       */
      pid_control_V1_DW.Memory1_PreviousInput[0] = pid_control_V1_B.Sum1[0];

      /* Update for Memory: '<S16>/Memory' incorporates:
       *  Sum: '<S16>/Sum'
       */
      pid_control_V1_DW.Memory_PreviousInput[1] = pid_control_V1_B.Sum[1];

      /* Update for Memory: '<S16>/Memory1' incorporates:
       *  Sum: '<S16>/Sum1'
       */
      pid_control_V1_DW.Memory1_PreviousInput[1] = pid_control_V1_B.Sum1[1];

      /* Update for Memory: '<S16>/Memory' incorporates:
       *  Sum: '<S16>/Sum'
       */
      pid_control_V1_DW.Memory_PreviousInput[2] = pid_control_V1_B.Sum[2];

      /* Update for Memory: '<S16>/Memory1' incorporates:
       *  Sum: '<S16>/Sum1'
       */
      pid_control_V1_DW.Memory1_PreviousInput[2] = pid_control_V1_B.Sum1[2];

      /* Update for RandomNumber: '<S288>/White Noise' */
      pid_control_V1_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
        (&pid_control_V1_DW.RandSeed);

      /* Update for RandomNumber: '<S306>/White Noise' */
      pid_control_V1_DW.NextOutput_j[0] = rt_nrand_Upu32_Yd_f_pw_snf
        (&pid_control_V1_DW.RandSeed_i[0]);
      pid_control_V1_DW.NextOutput_j[1] = rt_nrand_Upu32_Yd_f_pw_snf
        (&pid_control_V1_DW.RandSeed_i[1]);
      pid_control_V1_DW.NextOutput_j[2] = rt_nrand_Upu32_Yd_f_pw_snf
        (&pid_control_V1_DW.RandSeed_i[2]);
      pid_control_V1_DW.NextOutput_j[3] = rt_nrand_Upu32_Yd_f_pw_snf
        (&pid_control_V1_DW.RandSeed_i[3]);
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
  real_T tmp[2];
  _rtXdot = ((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs);

  /* Derivatives for Integrator: '<S16>/Integrator' */
  memcpy(&_rtXdot->Integrator_CSTATE[0], &pid_control_V1_B.XDOT[0], 12U * sizeof
         (real_T));

  /* Derivatives for Integrator: '<S107>/Integrator' */
  _rtXdot->Integrator_CSTATE_n = pid_control_V1_B.Switch_k;

  /* Derivatives for Integrator: '<S102>/Filter' */
  _rtXdot->Filter_CSTATE = pid_control_V1_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S53>/Integrator' */
  _rtXdot->Integrator_CSTATE_m = pid_control_V1_B.Switch_n;

  /* Derivatives for Integrator: '<S48>/Filter' */
  _rtXdot->Filter_CSTATE_g = pid_control_V1_B.FilterCoefficient_c;

  /* Derivatives for Integrator: '<S159>/Integrator' */
  _rtXdot->Integrator_CSTATE_p = pid_control_V1_B.SumI4;

  /* Derivatives for Integrator: '<S154>/Filter' */
  _rtXdot->Filter_CSTATE_m = pid_control_V1_B.FilterCoefficient_m;

  /* Derivatives for Integrator: '<S211>/Integrator' */
  _rtXdot->Integrator_CSTATE_d = pid_control_V1_B.IntegralGain;

  /* Derivatives for Integrator: '<S206>/Filter' */
  _rtXdot->Filter_CSTATE_f = pid_control_V1_B.FilterCoefficient_p;

  /* Derivatives for Integrator: '<S265>/Integrator' */
  _rtXdot->Integrator_CSTATE_f = pid_control_V1_B.Switch_j;

  /* Derivatives for Integrator: '<S260>/Filter' */
  _rtXdot->Filter_CSTATE_l = pid_control_V1_B.FilterCoefficient_cv;

  /* Derivatives for Integrator: '<S16>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = pid_control_V1_B.Power;

  /* Derivatives for Enabled SubSystem: '<S296>/Hpgw' */
  if (pid_control_V1_DW.Hpgw_MODE) {
    /* Derivatives for Integrator: '<S307>/pgw_p' */
    _rtXdot->pgw_p_CSTATE[0] = pid_control_V1_B.w_o[0];
    _rtXdot->pgw_p_CSTATE[1] = pid_control_V1_B.w_o[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->pgw_p_CSTATE[0]);
      for (i1=0; i1 < 2; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S296>/Hpgw' */

  /* Derivatives for Enabled SubSystem: '<S297>/Hwgw(s)' */
  if (pid_control_V1_DW.Hwgws_MODE) {
    /* Derivatives for Integrator: '<S312>/wg_p1' */
    _rtXdot->wg_p1_CSTATE[0] = pid_control_V1_B.w[0];

    /* Derivatives for Integrator: '<S312>/wg_p2' */
    _rtXdot->wg_p2_CSTATE[0] = pid_control_V1_B.w_a[0];

    /* Derivatives for Integrator: '<S312>/wg_p1' */
    _rtXdot->wg_p1_CSTATE[1] = pid_control_V1_B.w[1];

    /* Derivatives for Integrator: '<S312>/wg_p2' */
    _rtXdot->wg_p2_CSTATE[1] = pid_control_V1_B.w_a[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->wg_p1_CSTATE[0]);
      for (i1=0; i1 < 4; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S297>/Hwgw(s)' */

  /* Derivatives for Enabled SubSystem: '<S296>/Hqgw' */
  if (pid_control_V1_DW.Hqgw_MODE) {
    /* Derivatives for Integrator: '<S308>/qgw_p' */
    _rtXdot->qgw_p_CSTATE[0] = pid_control_V1_B.w_e0[0];
    _rtXdot->qgw_p_CSTATE[1] = pid_control_V1_B.w_e0[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->qgw_p_CSTATE[0]);
      for (i1=0; i1 < 2; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S296>/Hqgw' */

  /* Derivatives for Enabled SubSystem: '<S297>/Hvgw(s)' */
  if (pid_control_V1_DW.Hvgws_MODE) {
    /* Derivatives for Integrator: '<S311>/vg_p1' */
    _rtXdot->vg_p1_CSTATE[0] = pid_control_V1_B.w_g[0];

    /* Derivatives for Integrator: '<S311>/vgw_p2' */
    _rtXdot->vgw_p2_CSTATE[0] = pid_control_V1_B.w_e[0];

    /* Derivatives for Integrator: '<S311>/vg_p1' */
    _rtXdot->vg_p1_CSTATE[1] = pid_control_V1_B.w_g[1];

    /* Derivatives for Integrator: '<S311>/vgw_p2' */
    _rtXdot->vgw_p2_CSTATE[1] = pid_control_V1_B.w_e[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->vg_p1_CSTATE[0]);
      for (i1=0; i1 < 4; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S297>/Hvgw(s)' */

  /* Derivatives for Enabled SubSystem: '<S296>/Hrgw' */
  if (pid_control_V1_DW.Hrgw_MODE) {
    /* Derivatives for Integrator: '<S309>/rgw_p' */
    _rtXdot->rgw_p_CSTATE[0] = pid_control_V1_B.w_d[0];
    _rtXdot->rgw_p_CSTATE[1] = pid_control_V1_B.w_d[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->rgw_p_CSTATE[0]);
      for (i1=0; i1 < 2; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S296>/Hrgw' */

  /* Derivatives for Enabled SubSystem: '<S297>/Hugw(s)' */
  if (pid_control_V1_DW.Hugws_MODE) {
    /* Derivatives for Integrator: '<S310>/ug_p' */
    _rtXdot->ug_p_CSTATE[0] = pid_control_V1_B.w_n[0];
    _rtXdot->ug_p_CSTATE[1] = pid_control_V1_B.w_n[1];
  } else {
    {
      real_T *dx;
      int_T i1;
      dx = &(((XDot_pid_control_V1_T *) (&pid_control_V1_M)->derivs)
             ->ug_p_CSTATE[0]);
      for (i1=0; i1 < 2; i1++) {
        dx[i1] = 0.0;
      }
    }
  }

  /* End of Derivatives for SubSystem: '<S297>/Hugw(s)' */

  /* Derivatives for TransferFcn: '<S16>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE[0] = 0.0;
  _rtXdot->TransferFcn_CSTATE[0] += -0.898 *
    pid_control_V1_X.TransferFcn_CSTATE[0];
  _rtXdot->TransferFcn_CSTATE[1] = 0.0;
  _rtXdot->TransferFcn_CSTATE[0] += -0.806 *
    pid_control_V1_X.TransferFcn_CSTATE[1];
  _mm_storeu_pd(&tmp[0], _mm_add_pd(_mm_set_pd(_rtXdot->TransferFcn_CSTATE[0],
    pid_control_V1_X.TransferFcn_CSTATE[0]), _mm_set_pd(pid_control_V1_B.Output,
    _rtXdot->TransferFcn_CSTATE[1])));
  _rtXdot->TransferFcn_CSTATE[1] = tmp[0];
  _rtXdot->TransferFcn_CSTATE[0] = tmp[1];

  /* Derivatives for TransferFcn: '<S16>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += -0.01 * pid_control_V1_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += pid_control_V1_B.Output;
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

  /* Start for MATLABSystem: '<S13>/SourceBlock' */
  pid_control_V1_DW.obj_m.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_m.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_gq = true;
  pid_control_V1_DW.obj_m.isSetupComplete = false;
  pid_control_V1_DW.obj_m.isInitialized = 1;
  pid_con_Subscriber_setupImpl_on(&pid_control_V1_DW.obj_m);
  pid_control_V1_DW.obj_m.isSetupComplete = true;

  /* Start for MATLABSystem: '<S11>/SourceBlock' */
  pid_control_V1_DW.obj_f.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_f.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_f = true;
  pid_control_V1_DW.obj_f.isSetupComplete = false;
  pid_control_V1_DW.obj_f.isInitialized = 1;
  pid_contro_Subscriber_setupImpl(&pid_control_V1_DW.obj_f);
  pid_control_V1_DW.obj_f.isSetupComplete = true;

  /* Start for MATLABSystem: '<S14>/SourceBlock' */
  pid_control_V1_DW.obj_c.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_c.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_g = true;
  pid_control_V1_DW.obj_c.isSetupComplete = false;
  pid_control_V1_DW.obj_c.isInitialized = 1;
  pid_co_Subscriber_setupImpl_onh(&pid_control_V1_DW.obj_c);
  pid_control_V1_DW.obj_c.isSetupComplete = true;

  /* Start for MATLABSystem: '<S15>/SourceBlock' */
  pid_control_V1_DW.obj_n.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_n.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_k = true;
  pid_control_V1_DW.obj_n.isSetupComplete = false;
  pid_control_V1_DW.obj_n.isInitialized = 1;
  pid_c_Subscriber_setupImpl_onhg(&pid_control_V1_DW.obj_n);
  pid_control_V1_DW.obj_n.isSetupComplete = true;

  /* Start for MATLABSystem: '<S12>/SourceBlock' */
  pid_control_V1_DW.obj_k.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_k.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_g5 = true;
  pid_control_V1_DW.obj_k.isSetupComplete = false;
  pid_control_V1_DW.obj_k.isInitialized = 1;
  pid_cont_Subscriber_setupImpl_o(&pid_control_V1_DW.obj_k);
  pid_control_V1_DW.obj_k.isSetupComplete = true;

  /* Start for MATLABSystem: '<Root>/Coordinate Transformation Conversion' */
  pid_control_V1_DW.objisempty_d = true;
  pid_control_V1_DW.obj_cl.isInitialized = 1;

  /* Start for Atomic SubSystem: '<Root>/Call Service' */
  /* Start for MATLABSystem: '<S2>/ServiceCaller' */
  pid_control_V1_DW.obj.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_ft = true;
  pid_control_V1_DW.obj.isSetupComplete = false;
  pid_control_V1_DW.obj.isInitialized = 1;
  pid_con_ServiceCaller_setupImpl(&pid_control_V1_DW.obj);
  pid_control_V1_DW.obj.isSetupComplete = true;

  /* End of Start for SubSystem: '<Root>/Call Service' */

  /* Start for Enabled SubSystem: '<S296>/Hpgw' */
  (void) memset(&(pid_control_V1_XDis.pgw_p_CSTATE), 1,
                2*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S296>/Hpgw' */

  /* Start for Enabled SubSystem: '<S297>/Hwgw(s)' */
  (void) memset(&(pid_control_V1_XDis.wg_p1_CSTATE), 1,
                4*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S297>/Hwgw(s)' */

  /* Start for Enabled SubSystem: '<S296>/Hqgw' */
  (void) memset(&(pid_control_V1_XDis.qgw_p_CSTATE), 1,
                2*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S296>/Hqgw' */

  /* Start for Enabled SubSystem: '<S297>/Hvgw(s)' */
  (void) memset(&(pid_control_V1_XDis.vg_p1_CSTATE), 1,
                4*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S297>/Hvgw(s)' */

  /* Start for Enabled SubSystem: '<S296>/Hrgw' */
  (void) memset(&(pid_control_V1_XDis.rgw_p_CSTATE), 1,
                2*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S296>/Hrgw' */

  /* Start for Enabled SubSystem: '<S297>/Hugw(s)' */
  (void) memset(&(pid_control_V1_XDis.ug_p_CSTATE), 1,
                2*sizeof(boolean_T));

  /* End of Start for SubSystem: '<S297>/Hugw(s)' */

  /* Start for If: '<S301>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  pid_control_V1_DW.ifHeightMaxlowaltitudeelseifHei = -1;

  /* Start for If: '<S302>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  pid_control_V1_DW.ifHeightMaxlowaltitudeelseifH_k = -1;

  /* Start for MATLABSystem: '<S291>/SourceBlock' */
  pid_control_V1_DW.obj_hq.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_hq.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_a = true;
  pid_control_V1_DW.obj_hq.isSetupComplete = false;
  pid_control_V1_DW.obj_hq.isInitialized = 1;
  pid__Subscriber_setupImpl_onhgd(&pid_control_V1_DW.obj_hq);
  pid_control_V1_DW.obj_hq.isSetupComplete = true;

  /* Start for MATLABSystem: '<S292>/SourceBlock' */
  pid_control_V1_DW.obj_h4.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_h4.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_c = true;
  pid_control_V1_DW.obj_h4.isSetupComplete = false;
  pid_control_V1_DW.obj_h4.isInitialized = 1;
  pid_Subscriber_setupImpl_onhgd0(&pid_control_V1_DW.obj_h4);
  pid_control_V1_DW.obj_h4.isSetupComplete = true;

  /* Start for MATLABSystem: '<S293>/SourceBlock' */
  pid_control_V1_DW.obj_h.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_h.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty_l = true;
  pid_control_V1_DW.obj_h.isSetupComplete = false;
  pid_control_V1_DW.obj_h.isInitialized = 1;
  pi_Subscriber_setupImpl_onhgd03(&pid_control_V1_DW.obj_h);
  pid_control_V1_DW.obj_h.isSetupComplete = true;

  /* Start for MATLABSystem: '<S294>/SourceBlock' */
  pid_control_V1_DW.obj_p.QOSAvoidROSNamespaceConventions = false;
  pid_control_V1_DW.obj_p.matlabCodegenIsDeleted = false;
  pid_control_V1_DW.objisempty = true;
  pid_control_V1_DW.obj_p.isSetupComplete = false;
  pid_control_V1_DW.obj_p.isInitialized = 1;
  p_Subscriber_setupImpl_onhgd03r(&pid_control_V1_DW.obj_p);
  pid_control_V1_DW.obj_p.isSetupComplete = true;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay3' */
  pid_control_V1_DW.UnitDelay3_DSTATE = 1.0;

  /* InitializeConditions for Integrator: '<S16>/Integrator' */
  memcpy(&pid_control_V1_X.Integrator_CSTATE[0],
         &pid_control_V1_ConstP.Integrator_IC[0], 12U * sizeof(real_T));

  /* InitializeConditions for Integrator: '<S107>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_n = 0.0;

  /* InitializeConditions for Integrator: '<S102>/Filter' */
  pid_control_V1_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S53>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_m = 0.0;

  /* InitializeConditions for Integrator: '<S48>/Filter' */
  pid_control_V1_X.Filter_CSTATE_g = 0.0;

  /* InitializeConditions for Integrator: '<S159>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S154>/Filter' */
  pid_control_V1_X.Filter_CSTATE_m = 0.0;

  /* InitializeConditions for Integrator: '<S211>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_d = 0.0;

  /* InitializeConditions for Integrator: '<S206>/Filter' */
  pid_control_V1_X.Filter_CSTATE_f = 0.0;

  /* InitializeConditions for Integrator: '<S265>/Integrator' */
  pid_control_V1_X.Integrator_CSTATE_f = 0.0;

  /* InitializeConditions for Integrator: '<S260>/Filter' */
  pid_control_V1_X.Filter_CSTATE_l = 0.0;

  /* InitializeConditions for Integrator: '<S16>/Integrator1' */
  pid_control_V1_X.Integrator1_CSTATE = 0.0;

  /* InitializeConditions for RandomNumber: '<S288>/White Noise' */
  pid_control_V1_DW.RandSeed = 1529675776U;
  pid_control_V1_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
    (&pid_control_V1_DW.RandSeed);

  /* InitializeConditions for RandomNumber: '<S306>/White Noise' */
  pid_control_V1_DW.RandSeed_i[0] = 1529675776U;
  pid_control_V1_DW.NextOutput_j[0] = rt_nrand_Upu32_Yd_f_pw_snf
    (&pid_control_V1_DW.RandSeed_i[0]);
  pid_control_V1_DW.RandSeed_i[1] = 1529741312U;
  pid_control_V1_DW.NextOutput_j[1] = rt_nrand_Upu32_Yd_f_pw_snf
    (&pid_control_V1_DW.RandSeed_i[1]);
  pid_control_V1_DW.RandSeed_i[2] = 1529806848U;
  pid_control_V1_DW.NextOutput_j[2] = rt_nrand_Upu32_Yd_f_pw_snf
    (&pid_control_V1_DW.RandSeed_i[2]);
  pid_control_V1_DW.RandSeed_i[3] = 1529872384U;
  pid_control_V1_DW.NextOutput_j[3] = rt_nrand_Upu32_Yd_f_pw_snf
    (&pid_control_V1_DW.RandSeed_i[3]);

  /* InitializeConditions for TransferFcn: '<S16>/Transfer Fcn' */
  pid_control_V1_X.TransferFcn_CSTATE[0] = 0.0;
  pid_control_V1_X.TransferFcn_CSTATE[1] = 0.0;

  /* InitializeConditions for TransferFcn: '<S16>/Transfer Fcn1' */
  pid_control_V1_X.TransferFcn1_CSTATE = 0.0;

  /* SystemInitialize for Enabled SubSystem: '<S13>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_b);

  /* End of SystemInitialize for SubSystem: '<S13>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S11>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem);

  /* End of SystemInitialize for SubSystem: '<S11>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S14>/Enabled Subsystem' */
  pid_con_EnabledSubsystem_d_Init(&pid_control_V1_B.EnabledSubsystem_h);

  /* End of SystemInitialize for SubSystem: '<S14>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S15>/Enabled Subsystem' */
  pid_con_EnabledSubsystem_d_Init(&pid_control_V1_B.EnabledSubsystem_bk);

  /* End of SystemInitialize for SubSystem: '<S15>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S12>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_a);

  /* End of SystemInitialize for SubSystem: '<S12>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hpgw' */
  /* InitializeConditions for Integrator: '<S307>/pgw_p' */
  pid_control_V1_X.pgw_p_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hpgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hwgw(s)' */
  /* InitializeConditions for Integrator: '<S312>/wg_p1' */
  pid_control_V1_X.wg_p1_CSTATE[0] = 0.0;

  /* InitializeConditions for Integrator: '<S312>/wg_p2' */
  pid_control_V1_X.wg_p2_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hwgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hqgw' */
  /* InitializeConditions for Integrator: '<S308>/qgw_p' */
  pid_control_V1_X.qgw_p_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hqgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hvgw(s)' */
  /* InitializeConditions for Integrator: '<S311>/vg_p1' */
  pid_control_V1_X.vg_p1_CSTATE[0] = 0.0;

  /* InitializeConditions for Integrator: '<S311>/vgw_p2' */
  pid_control_V1_X.vgw_p2_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hvgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hrgw' */
  /* InitializeConditions for Integrator: '<S309>/rgw_p' */
  pid_control_V1_X.rgw_p_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hrgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hugw(s)' */
  /* InitializeConditions for Integrator: '<S310>/ug_p' */
  pid_control_V1_X.ug_p_CSTATE[0] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hugw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hpgw' */
  /* InitializeConditions for Integrator: '<S307>/pgw_p' */
  pid_control_V1_X.pgw_p_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hpgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hwgw(s)' */
  /* InitializeConditions for Integrator: '<S312>/wg_p1' */
  pid_control_V1_X.wg_p1_CSTATE[1] = 0.0;

  /* InitializeConditions for Integrator: '<S312>/wg_p2' */
  pid_control_V1_X.wg_p2_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hwgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hqgw' */
  /* InitializeConditions for Integrator: '<S308>/qgw_p' */
  pid_control_V1_X.qgw_p_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hqgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hvgw(s)' */
  /* InitializeConditions for Integrator: '<S311>/vg_p1' */
  pid_control_V1_X.vg_p1_CSTATE[1] = 0.0;

  /* InitializeConditions for Integrator: '<S311>/vgw_p2' */
  pid_control_V1_X.vgw_p2_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hvgw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S296>/Hrgw' */
  /* InitializeConditions for Integrator: '<S309>/rgw_p' */
  pid_control_V1_X.rgw_p_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S296>/Hrgw' */

  /* SystemInitialize for Enabled SubSystem: '<S297>/Hugw(s)' */
  /* InitializeConditions for Integrator: '<S310>/ug_p' */
  pid_control_V1_X.ug_p_CSTATE[1] = 0.0;

  /* End of SystemInitialize for SubSystem: '<S297>/Hugw(s)' */

  /* SystemInitialize for Enabled SubSystem: '<S291>/Enabled Subsystem' */
  pid_con_EnabledSubsystem_i_Init(&pid_control_V1_B.EnabledSubsystem_p);

  /* End of SystemInitialize for SubSystem: '<S291>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S292>/Enabled Subsystem' */
  pid_con_EnabledSubsystem_i_Init(&pid_control_V1_B.EnabledSubsystem_g);

  /* End of SystemInitialize for SubSystem: '<S292>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S293>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_k);

  /* End of SystemInitialize for SubSystem: '<S293>/Enabled Subsystem' */

  /* SystemInitialize for Enabled SubSystem: '<S294>/Enabled Subsystem' */
  pid_contr_EnabledSubsystem_Init(&pid_control_V1_B.EnabledSubsystem_pu);

  /* End of SystemInitialize for SubSystem: '<S294>/Enabled Subsystem' */
}

/* Model terminate function */
void pid_control_V1::terminate()
{
  /* Terminate for MATLABSystem: '<S13>/SourceBlock' */
  if (!pid_control_V1_DW.obj_m.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_m.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_m.isInitialized == 1) &&
        pid_control_V1_DW.obj_m.isSetupComplete) {
      Sub_pid_control_V1_435.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S13>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S11>/SourceBlock' */
  if (!pid_control_V1_DW.obj_f.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_f.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_f.isInitialized == 1) &&
        pid_control_V1_DW.obj_f.isSetupComplete) {
      Sub_pid_control_V1_466.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S11>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S14>/SourceBlock' */
  if (!pid_control_V1_DW.obj_c.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_c.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_c.isInitialized == 1) &&
        pid_control_V1_DW.obj_c.isSetupComplete) {
      Sub_pid_control_V1_467.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S14>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S15>/SourceBlock' */
  if (!pid_control_V1_DW.obj_n.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_n.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_n.isInitialized == 1) &&
        pid_control_V1_DW.obj_n.isSetupComplete) {
      Sub_pid_control_V1_476.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S15>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S12>/SourceBlock' */
  if (!pid_control_V1_DW.obj_k.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_k.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_k.isInitialized == 1) &&
        pid_control_V1_DW.obj_k.isSetupComplete) {
      Sub_pid_control_V1_377.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S12>/SourceBlock' */

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
  /* Terminate for MATLABSystem: '<S291>/SourceBlock' */
  if (!pid_control_V1_DW.obj_hq.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_hq.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_hq.isInitialized == 1) &&
        pid_control_V1_DW.obj_hq.isSetupComplete) {
      Sub_pid_control_V1_417.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S291>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S292>/SourceBlock' */
  if (!pid_control_V1_DW.obj_h4.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_h4.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_h4.isInitialized == 1) &&
        pid_control_V1_DW.obj_h4.isSetupComplete) {
      Sub_pid_control_V1_423.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S292>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S293>/SourceBlock' */
  if (!pid_control_V1_DW.obj_h.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_h.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_h.isInitialized == 1) &&
        pid_control_V1_DW.obj_h.isSetupComplete) {
      Sub_pid_control_V1_443.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S293>/SourceBlock' */

  /* Terminate for MATLABSystem: '<S294>/SourceBlock' */
  if (!pid_control_V1_DW.obj_p.matlabCodegenIsDeleted) {
    pid_control_V1_DW.obj_p.matlabCodegenIsDeleted = true;
    if ((pid_control_V1_DW.obj_p.isInitialized == 1) &&
        pid_control_V1_DW.obj_p.isSetupComplete) {
      Sub_pid_control_V1_445.resetSubscriberPtr();//();
    }
  }

  /* End of Terminate for MATLABSystem: '<S294>/SourceBlock' */
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
