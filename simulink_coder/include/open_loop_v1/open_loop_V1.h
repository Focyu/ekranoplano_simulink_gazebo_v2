/*
 * open_loop_V1.h
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

#ifndef open_loop_V1_h_
#define open_loop_V1_h_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "slros2_initialize.h"
#include "open_loop_V1_types.h"

extern "C"
{

#include "rt_nonfinite.h"

}

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

#include <string.h>
#include <stddef.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

#ifndef rtmGetTStart
#define rtmGetTStart(rtm)              ((rtm)->Timing.tStart)
#endif

/* Block signals (default storage) */
struct B_open_loop_V1_T {
  SL_Bus_gazebo_msgs_SetEntityStateRequest BusAssignment;/* '<Root>/Bus Assignment' */
  real_T x[12];                        /* '<S7>/Integrator' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function1' */
  uint8_T stringOut_l[128];            /* '<Root>/MATLAB Function' */
  real_T c_theta_g[9];
  real_T control_vector[5];
  char_T b_zeroDelimTopic[25];
  real_T Vb_w[3];
  real_T wbe_b[3];
  real_T F_b[3];
  real_T CD_iw_IGE_m[3];
  real_T Mcg_b[3];
  real_T dv[2];
  real_T FilterCoefficient;            /* '<S97>/Filter Coefficient' */
  real_T Saturation;                   /* '<S101>/Saturation' */
  real_T FilterCoefficient_m;          /* '<S45>/Filter Coefficient' */
  real_T Saturation_f;                 /* '<S49>/Saturation' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T Power;                        /* '<S7>/Product2' */
  real_T Gain3;                        /* '<S7>/Gain3' */
  real_T powerdemand;                  /* '<S7>/Divide' */
  real_T loadtorque;                   /* '<S7>/Divide1' */
  real_T IntegralGain;                 /* '<S39>/Integral Gain' */
  real_T IntegralGain_m;               /* '<S91>/Integral Gain' */
  real_T EnergykWh;                    /* '<S7>/Gain1' */
  real_T XDOT[40];                     /* '<S7>/MATLAB Function' */
  real_T u2;
  real_T u4;
  real_T u5;
  real_T c_phi;
  real_T s_phi;
  real_T c_theta;
  real_T s_theta;
  real_T c_psi;
  real_T s_psi;
  real_T Va;
  real_T CL_w_OGE;
  real_T CL_h_OGE;
  real_T CL_w_IGE;
  real_T CD_iw_IGE;
  real_T CD_ih_IGE;
  real_T Cl;
  real_T Cm;
  real_T Cn;
  real_T FA_b_tmp;
  real_T FA_b_tmp_m;
  real_T Va_b_idx_2;
  real_T Mcg_b_idx_1;
  real_T Mcg_b_idx_2;
  real_T FA_b_idx_0;
  real_T FA_b_idx_1;
  real_T FA_b_idx_2;
  real_T c_theta_tmp;
  real_T c_theta_tmp_c;
  real_T c_theta_tmp_k;
  real_T c_theta_tmp_cx;
  real_T c_theta_tmp_b;
  real_T c_theta_tmp_p;
  real_T c_theta_tmp_cv;
  real_T c_theta_tmp_f;
  real_T d;
  real_T d1;
  int32_T i;
  int32_T i_g;
  int32_T i1;
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function1' */
  uint32_T lengthOut_e;                /* '<Root>/MATLAB Function' */
  boolean_T serverAvailableOnTime;
  boolean_T b;
  SL_Bus_gazebo_msgs_SetEntityStateResponse r;
};

/* Block states (default storage) for system '<Root>' */
struct DW_open_loop_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  struct {
    void *LoggedData;
  } ToWorkspace10_PWORK;               /* '<S7>/To Workspace10' */

  struct {
    void *LoggedData;
  } ToWorkspace11_PWORK;               /* '<S7>/To Workspace11' */

  struct {
    void *LoggedData;
  } ToWorkspace12_PWORK;               /* '<S7>/To Workspace12' */

  struct {
    void *LoggedData;
  } ToWorkspace13_PWORK;               /* '<S7>/To Workspace13' */

  struct {
    void *LoggedData;
  } ToWorkspace14_PWORK;               /* '<S7>/To Workspace14' */

  struct {
    void *LoggedData;
  } ToWorkspace15_PWORK;               /* '<S7>/To Workspace15' */

  struct {
    void *LoggedData;
  } ToWorkspace16_PWORK;               /* '<S7>/To Workspace16' */

  struct {
    void *LoggedData;
  } ToWorkspace17_PWORK;               /* '<S7>/To Workspace17' */

  struct {
    void *LoggedData;
  } ToWorkspace18_PWORK;               /* '<S7>/To Workspace18' */

  struct {
    void *LoggedData;
  } ToWorkspace19_PWORK;               /* '<S7>/To Workspace19' */

  struct {
    void *LoggedData;
  } ToWorkspace2_PWORK;                /* '<S7>/To Workspace2' */

  struct {
    void *LoggedData;
  } ToWorkspace20_PWORK;               /* '<S7>/To Workspace20' */

  struct {
    void *LoggedData;
  } ToWorkspace21_PWORK;               /* '<S7>/To Workspace21' */

  struct {
    void *LoggedData;
  } ToWorkspace22_PWORK;               /* '<S7>/To Workspace22' */

  struct {
    void *LoggedData;
  } ToWorkspace23_PWORK;               /* '<S7>/To Workspace23' */

  struct {
    void *LoggedData;
  } ToWorkspace24_PWORK;               /* '<S7>/To Workspace24' */

  struct {
    void *LoggedData;
  } ToWorkspace25_PWORK;               /* '<S7>/To Workspace25' */

  struct {
    void *LoggedData;
  } ToWorkspace26_PWORK;               /* '<S7>/To Workspace26' */

  struct {
    void *LoggedData;
  } ToWorkspace27_PWORK;               /* '<S7>/To Workspace27' */

  struct {
    void *LoggedData;
  } ToWorkspace28_PWORK;               /* '<S7>/To Workspace28' */

  struct {
    void *LoggedData;
  } ToWorkspace29_PWORK;               /* '<S7>/To Workspace29' */

  struct {
    void *LoggedData;
  } ToWorkspace3_PWORK;                /* '<S7>/To Workspace3' */

  struct {
    void *LoggedData;
  } ToWorkspace4_PWORK;                /* '<S7>/To Workspace4' */

  struct {
    void *LoggedData;
  } ToWorkspace5_PWORK;                /* '<S7>/To Workspace5' */

  struct {
    void *LoggedData;
  } ToWorkspace6_PWORK;                /* '<S7>/To Workspace6' */

  struct {
    void *LoggedData;
  } ToWorkspace7_PWORK;                /* '<S7>/To Workspace7' */

  struct {
    void *LoggedData;
  } ToWorkspace8_PWORK;                /* '<S7>/To Workspace8' */

  struct {
    void *LoggedData;
  } ToWorkspace9_PWORK;                /* '<S7>/To Workspace9' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<S7>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S7>/To Workspace1' */

  struct {
    void *LoggedData;
  } ToWorkspace30_PWORK;               /* '<S7>/To Workspace30' */

  struct {
    void *LoggedData;
  } ToWorkspace32_PWORK;               /* '<S7>/To Workspace32' */

  struct {
    void *LoggedData;
  } ToWorkspace33_PWORK;               /* '<S7>/To Workspace33' */

  struct {
    void *LoggedData;
  } ToWorkspace31_PWORK;               /* '<S7>/To Workspace31' */

  robotics_slcore_internal_bloc_T obj_c;
                             /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty;      /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_f;              /* '<S2>/ServiceCaller' */
};

/* Continuous states (default storage) */
struct X_open_loop_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S7>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S94>/Integrator' */
  real_T Filter_CSTATE;                /* '<S89>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S42>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S37>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S7>/Integrator1' */
};

/* State derivatives (default storage) */
struct XDot_open_loop_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S7>/Integrator' */
  real_T Integrator_CSTATE_b;          /* '<S94>/Integrator' */
  real_T Filter_CSTATE;                /* '<S89>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S42>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S37>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S7>/Integrator1' */
};

/* State disabled  */
struct XDis_open_loop_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S7>/Integrator' */
  boolean_T Integrator_CSTATE_b;       /* '<S94>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S89>/Filter' */
  boolean_T Integrator_CSTATE_p;       /* '<S42>/Integrator' */
  boolean_T Filter_CSTATE_m;           /* '<S37>/Filter' */
  boolean_T Integrator1_CSTATE;        /* '<S7>/Integrator1' */
};

/* Invariant block signals (default storage) */
struct ConstB_open_loop_V1_T {
  real_T motorspeed;                   /* '<S7>/Gain2' */
};

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
struct ODE3_IntgData {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
};

#endif

/* Constant parameters (default storage) */
struct ConstP_open_loop_V1_T {
  /* Expression: x0
   * Referenced by: '<S7>/Integrator'
   */
  real_T Integrator_IC[12];
};

/* Real-time Model Data Structure */
struct tag_RTM_open_loop_V1_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_open_loop_V1_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_open_loop_V1_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[17];
  real_T odeF[3][17];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tStart;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

extern const ConstB_open_loop_V1_T open_loop_V1_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_open_loop_V1_T open_loop_V1_ConstP;

/* Class declaration for model open_loop_V1 */
class open_loop_V1
{
  /* public data and function members */
 public:
  /* Real-Time Model get method */
  RT_MODEL_open_loop_V1_T * getRTM();

  /* model start function */
  void start();

  /* Initial conditions function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  void terminate();

  /* Constructor */
  open_loop_V1();

  /* Destructor */
  ~open_loop_V1();

  /* private data and function members */
 private:
  /* Block signals */
  B_open_loop_V1_T open_loop_V1_B;

  /* Block states */
  DW_open_loop_V1_T open_loop_V1_DW;

  /* Block continuous states */
  X_open_loop_V1_T open_loop_V1_X;

  /* Block Continuous state disabled vector */
  XDis_open_loop_V1_T open_loop_V1_XDis;

  /* private member function(s) for subsystem '<Root>'*/
  void open_lo_ServiceCaller_setupImpl(const ros_slros2_internal_block_Ser_T
    *obj);
  real_T open_loop_V1_rt_atan2d_snf(real_T u0, real_T u1);
  real_T open_loop_V1_rt_powd_snf(real_T u0, real_T u1);

  /* Global mass matrix */

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void open_loop_V1_derivatives();

  /* Real-Time Model */
  RT_MODEL_open_loop_V1_T open_loop_V1_M;
};

extern volatile boolean_T stopRequested;
extern volatile boolean_T runModel;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/Display' : Unused code path elimination
 * Block '<Root>/Gain4' : Eliminated nontunable gain of 1
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'open_loop_V1'
 * '<S1>'   : 'open_loop_V1/Blank Message'
 * '<S2>'   : 'open_loop_V1/Call Service'
 * '<S3>'   : 'open_loop_V1/MATLAB Function'
 * '<S4>'   : 'open_loop_V1/MATLAB Function1'
 * '<S5>'   : 'open_loop_V1/PID Controller'
 * '<S6>'   : 'open_loop_V1/PID Controller1'
 * '<S7>'   : 'open_loop_V1/Subsystem'
 * '<S8>'   : 'open_loop_V1/PID Controller/Anti-windup'
 * '<S9>'   : 'open_loop_V1/PID Controller/D Gain'
 * '<S10>'  : 'open_loop_V1/PID Controller/External Derivative'
 * '<S11>'  : 'open_loop_V1/PID Controller/Filter'
 * '<S12>'  : 'open_loop_V1/PID Controller/Filter ICs'
 * '<S13>'  : 'open_loop_V1/PID Controller/I Gain'
 * '<S14>'  : 'open_loop_V1/PID Controller/Ideal P Gain'
 * '<S15>'  : 'open_loop_V1/PID Controller/Ideal P Gain Fdbk'
 * '<S16>'  : 'open_loop_V1/PID Controller/Integrator'
 * '<S17>'  : 'open_loop_V1/PID Controller/Integrator ICs'
 * '<S18>'  : 'open_loop_V1/PID Controller/N Copy'
 * '<S19>'  : 'open_loop_V1/PID Controller/N Gain'
 * '<S20>'  : 'open_loop_V1/PID Controller/P Copy'
 * '<S21>'  : 'open_loop_V1/PID Controller/Parallel P Gain'
 * '<S22>'  : 'open_loop_V1/PID Controller/Reset Signal'
 * '<S23>'  : 'open_loop_V1/PID Controller/Saturation'
 * '<S24>'  : 'open_loop_V1/PID Controller/Saturation Fdbk'
 * '<S25>'  : 'open_loop_V1/PID Controller/Sum'
 * '<S26>'  : 'open_loop_V1/PID Controller/Sum Fdbk'
 * '<S27>'  : 'open_loop_V1/PID Controller/Tracking Mode'
 * '<S28>'  : 'open_loop_V1/PID Controller/Tracking Mode Sum'
 * '<S29>'  : 'open_loop_V1/PID Controller/Tsamp - Integral'
 * '<S30>'  : 'open_loop_V1/PID Controller/Tsamp - Ngain'
 * '<S31>'  : 'open_loop_V1/PID Controller/postSat Signal'
 * '<S32>'  : 'open_loop_V1/PID Controller/preInt Signal'
 * '<S33>'  : 'open_loop_V1/PID Controller/preSat Signal'
 * '<S34>'  : 'open_loop_V1/PID Controller/Anti-windup/Passthrough'
 * '<S35>'  : 'open_loop_V1/PID Controller/D Gain/Internal Parameters'
 * '<S36>'  : 'open_loop_V1/PID Controller/External Derivative/Error'
 * '<S37>'  : 'open_loop_V1/PID Controller/Filter/Cont. Filter'
 * '<S38>'  : 'open_loop_V1/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S39>'  : 'open_loop_V1/PID Controller/I Gain/Internal Parameters'
 * '<S40>'  : 'open_loop_V1/PID Controller/Ideal P Gain/Passthrough'
 * '<S41>'  : 'open_loop_V1/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S42>'  : 'open_loop_V1/PID Controller/Integrator/Continuous'
 * '<S43>'  : 'open_loop_V1/PID Controller/Integrator ICs/Internal IC'
 * '<S44>'  : 'open_loop_V1/PID Controller/N Copy/Disabled'
 * '<S45>'  : 'open_loop_V1/PID Controller/N Gain/Internal Parameters'
 * '<S46>'  : 'open_loop_V1/PID Controller/P Copy/Disabled'
 * '<S47>'  : 'open_loop_V1/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S48>'  : 'open_loop_V1/PID Controller/Reset Signal/Disabled'
 * '<S49>'  : 'open_loop_V1/PID Controller/Saturation/Enabled'
 * '<S50>'  : 'open_loop_V1/PID Controller/Saturation Fdbk/Disabled'
 * '<S51>'  : 'open_loop_V1/PID Controller/Sum/Sum_PID'
 * '<S52>'  : 'open_loop_V1/PID Controller/Sum Fdbk/Disabled'
 * '<S53>'  : 'open_loop_V1/PID Controller/Tracking Mode/Disabled'
 * '<S54>'  : 'open_loop_V1/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S55>'  : 'open_loop_V1/PID Controller/Tsamp - Integral/TsSignalSpecification'
 * '<S56>'  : 'open_loop_V1/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S57>'  : 'open_loop_V1/PID Controller/postSat Signal/Forward_Path'
 * '<S58>'  : 'open_loop_V1/PID Controller/preInt Signal/Internal PreInt'
 * '<S59>'  : 'open_loop_V1/PID Controller/preSat Signal/Forward_Path'
 * '<S60>'  : 'open_loop_V1/PID Controller1/Anti-windup'
 * '<S61>'  : 'open_loop_V1/PID Controller1/D Gain'
 * '<S62>'  : 'open_loop_V1/PID Controller1/External Derivative'
 * '<S63>'  : 'open_loop_V1/PID Controller1/Filter'
 * '<S64>'  : 'open_loop_V1/PID Controller1/Filter ICs'
 * '<S65>'  : 'open_loop_V1/PID Controller1/I Gain'
 * '<S66>'  : 'open_loop_V1/PID Controller1/Ideal P Gain'
 * '<S67>'  : 'open_loop_V1/PID Controller1/Ideal P Gain Fdbk'
 * '<S68>'  : 'open_loop_V1/PID Controller1/Integrator'
 * '<S69>'  : 'open_loop_V1/PID Controller1/Integrator ICs'
 * '<S70>'  : 'open_loop_V1/PID Controller1/N Copy'
 * '<S71>'  : 'open_loop_V1/PID Controller1/N Gain'
 * '<S72>'  : 'open_loop_V1/PID Controller1/P Copy'
 * '<S73>'  : 'open_loop_V1/PID Controller1/Parallel P Gain'
 * '<S74>'  : 'open_loop_V1/PID Controller1/Reset Signal'
 * '<S75>'  : 'open_loop_V1/PID Controller1/Saturation'
 * '<S76>'  : 'open_loop_V1/PID Controller1/Saturation Fdbk'
 * '<S77>'  : 'open_loop_V1/PID Controller1/Sum'
 * '<S78>'  : 'open_loop_V1/PID Controller1/Sum Fdbk'
 * '<S79>'  : 'open_loop_V1/PID Controller1/Tracking Mode'
 * '<S80>'  : 'open_loop_V1/PID Controller1/Tracking Mode Sum'
 * '<S81>'  : 'open_loop_V1/PID Controller1/Tsamp - Integral'
 * '<S82>'  : 'open_loop_V1/PID Controller1/Tsamp - Ngain'
 * '<S83>'  : 'open_loop_V1/PID Controller1/postSat Signal'
 * '<S84>'  : 'open_loop_V1/PID Controller1/preInt Signal'
 * '<S85>'  : 'open_loop_V1/PID Controller1/preSat Signal'
 * '<S86>'  : 'open_loop_V1/PID Controller1/Anti-windup/Passthrough'
 * '<S87>'  : 'open_loop_V1/PID Controller1/D Gain/Internal Parameters'
 * '<S88>'  : 'open_loop_V1/PID Controller1/External Derivative/Error'
 * '<S89>'  : 'open_loop_V1/PID Controller1/Filter/Cont. Filter'
 * '<S90>'  : 'open_loop_V1/PID Controller1/Filter ICs/Internal IC - Filter'
 * '<S91>'  : 'open_loop_V1/PID Controller1/I Gain/Internal Parameters'
 * '<S92>'  : 'open_loop_V1/PID Controller1/Ideal P Gain/Passthrough'
 * '<S93>'  : 'open_loop_V1/PID Controller1/Ideal P Gain Fdbk/Disabled'
 * '<S94>'  : 'open_loop_V1/PID Controller1/Integrator/Continuous'
 * '<S95>'  : 'open_loop_V1/PID Controller1/Integrator ICs/Internal IC'
 * '<S96>'  : 'open_loop_V1/PID Controller1/N Copy/Disabled'
 * '<S97>'  : 'open_loop_V1/PID Controller1/N Gain/Internal Parameters'
 * '<S98>'  : 'open_loop_V1/PID Controller1/P Copy/Disabled'
 * '<S99>'  : 'open_loop_V1/PID Controller1/Parallel P Gain/Internal Parameters'
 * '<S100>' : 'open_loop_V1/PID Controller1/Reset Signal/Disabled'
 * '<S101>' : 'open_loop_V1/PID Controller1/Saturation/Enabled'
 * '<S102>' : 'open_loop_V1/PID Controller1/Saturation Fdbk/Disabled'
 * '<S103>' : 'open_loop_V1/PID Controller1/Sum/Sum_PID'
 * '<S104>' : 'open_loop_V1/PID Controller1/Sum Fdbk/Disabled'
 * '<S105>' : 'open_loop_V1/PID Controller1/Tracking Mode/Disabled'
 * '<S106>' : 'open_loop_V1/PID Controller1/Tracking Mode Sum/Passthrough'
 * '<S107>' : 'open_loop_V1/PID Controller1/Tsamp - Integral/TsSignalSpecification'
 * '<S108>' : 'open_loop_V1/PID Controller1/Tsamp - Ngain/Passthrough'
 * '<S109>' : 'open_loop_V1/PID Controller1/postSat Signal/Forward_Path'
 * '<S110>' : 'open_loop_V1/PID Controller1/preInt Signal/Internal PreInt'
 * '<S111>' : 'open_loop_V1/PID Controller1/preSat Signal/Forward_Path'
 * '<S112>' : 'open_loop_V1/Subsystem/MATLAB Function'
 */
#endif                                 /* open_loop_V1_h_ */
