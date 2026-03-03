/*
 * pid_control_V1.h
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

#ifndef pid_control_V1_h_
#define pid_control_V1_h_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "slros2_initialize.h"
#include "pid_control_V1_types.h"

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

#include <string.h>

extern "C"
{

#include "rt_nonfinite.h"

}

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

/* Block signals for system '<S9>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_contro_T {
  SL_Bus_std_msgs_Float64 In1;         /* '<S225>/In1' */
};

/* Block signals (default storage) */
struct B_pid_control_V1_T {
  SL_Bus_gazebo_msgs_SetEntityStateRequest BusAssignment;/* '<Root>/Bus Assignment' */
  real_T x[12];                        /* '<S12>/Integrator' */
  real_T FA_b_tmp[9];
  real_T TmpSignalConversionAtSFunct[5];/* '<S12>/MATLAB Function' */
  char_T b_zeroDelimTopic[25];
  real_T wbe_b[3];
  real_T FE2_b[3];
  real_T F_b[3];
  real_T Q[3];
  char_T b_zeroDelimTopic_m[17];
  char_T b_zeroDelimTopic_c[16];
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  real_T dv[2];
  real_T Switch1;                      /* '<Root>/Switch1' */
  real_T FilterCoefficient;            /* '<S104>/Filter Coefficient' */
  real_T Saturation;                   /* '<S108>/Saturation' */
  real_T Switch2;                      /* '<Root>/Switch2' */
  real_T FilterCoefficient_p;          /* '<S210>/Filter Coefficient' */
  real_T Saturation_m;                 /* '<S214>/Saturation' */
  real_T FilterCoefficient_c;          /* '<S50>/Filter Coefficient' */
  real_T Saturation_k;                 /* '<S54>/Saturation' */
  real_T Switch;                       /* '<Root>/Switch' */
  real_T FilterCoefficient_g;          /* '<S158>/Filter Coefficient' */
  real_T Saturation_p;                 /* '<S162>/Saturation' */
  real_T Gain5;                        /* '<Root>/Gain5' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T Power;                        /* '<S12>/Product2' */
  real_T Gain3;                        /* '<S12>/Gain3' */
  real_T powerdemand;                  /* '<S12>/Divide' */
  real_T loadtorque;                   /* '<S12>/Divide1' */
  real_T Switch_c;                     /* '<S145>/Switch' */
  real_T IntegralGain;                 /* '<S44>/Integral Gain' */
  real_T IntegralGain_a;               /* '<S204>/Integral Gain' */
  real_T Switch_a;                     /* '<S91>/Switch' */
  real_T EnergykWh;                    /* '<S12>/Gain1' */
  real_T XDOT[40];                     /* '<S12>/MATLAB Function' */
  real_T u2;
  real_T Va;
  real_T alpha;
  real_T beta;
  real_T hw;
  real_T Q_k;
  real_T CL_w_OGE;
  real_T CL_h_OGE;
  real_T CL_h_IGE;
  real_T CD_iw_IGE;
  real_T CD_ih_IGE;
  real_T CQ;
  real_T Cl;
  real_T Cn;
  real_T Vd1;
  real_T Tp1;
  real_T SignPreSat;                   /* '<S91>/SignPreSat' */
  real_T Sum5;                         /* '<Root>/Sum5' */
  real_T SignPreSat_o;                 /* '<S145>/SignPreSat' */
  real_T FE_b;
  real_T FE1_b_idx_0;
  real_T FE1_b_idx_1;
  real_T FE1_b_idx_2;
  real_T Mcg_b_idx_0;
  real_T Mcg_b_idx_2;
  real_T FE_b_idx_0;
  real_T FA_b_idx_0;
  real_T FA_b_idx_1;
  real_T FA_b_idx_2;
  real_T FE2_b_c;
  real_T FE2_b_b;
  real_T FE2_b_p;
  SL_Bus_std_msgs_Float64 SourceBlock_o2_l;/* '<S9>/SourceBlock' */
  SL_Bus_std_msgs_Float64 SourceBlock_o2;/* '<S11>/SourceBlock' */
  SL_Bus_std_msgs_Float64 SourceBlock_o2_e;/* '<S10>/SourceBlock' */
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function1' */
  uint32_T lengthOut_e;                /* '<Root>/MATLAB Function' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function1' */
  uint8_T stringOut_l[128];            /* '<Root>/MATLAB Function' */
  boolean_T AND3;                      /* '<S145>/AND3' */
  boolean_T Memory;                    /* '<S145>/Memory' */
  boolean_T AND3_j;                    /* '<S91>/AND3' */
  boolean_T Memory_n;                  /* '<S91>/Memory' */
  boolean_T SourceBlock_o1;            /* '<S11>/SourceBlock' */
  boolean_T SourceBlock_o1_l;          /* '<S10>/SourceBlock' */
  boolean_T SourceBlock_o1_i;          /* '<S9>/SourceBlock' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_a;/* '<S11>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_c;/* '<S10>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem;/* '<S9>/Enabled Subsystem' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_pid_control_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  ros_slros2_internal_block_Sub_T obj_k;/* '<S11>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_j;/* '<S10>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_h;/* '<S9>/SourceBlock' */
  real_T UnitDelay1_DSTATE;            /* '<Root>/Unit Delay1' */
  real_T UnitDelay2_DSTATE;            /* '<Root>/Unit Delay2' */
  real_T UnitDelay_DSTATE;             /* '<Root>/Unit Delay' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<Root>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace10_PWORK;               /* '<S12>/To Workspace10' */

  struct {
    void *LoggedData;
  } ToWorkspace11_PWORK;               /* '<S12>/To Workspace11' */

  struct {
    void *LoggedData;
  } ToWorkspace12_PWORK;               /* '<S12>/To Workspace12' */

  struct {
    void *LoggedData;
  } ToWorkspace13_PWORK;               /* '<S12>/To Workspace13' */

  struct {
    void *LoggedData;
  } ToWorkspace14_PWORK;               /* '<S12>/To Workspace14' */

  struct {
    void *LoggedData;
  } ToWorkspace15_PWORK;               /* '<S12>/To Workspace15' */

  struct {
    void *LoggedData;
  } ToWorkspace16_PWORK;               /* '<S12>/To Workspace16' */

  struct {
    void *LoggedData;
  } ToWorkspace17_PWORK;               /* '<S12>/To Workspace17' */

  struct {
    void *LoggedData;
  } ToWorkspace18_PWORK;               /* '<S12>/To Workspace18' */

  struct {
    void *LoggedData;
  } ToWorkspace19_PWORK;               /* '<S12>/To Workspace19' */

  struct {
    void *LoggedData;
  } ToWorkspace2_PWORK;                /* '<S12>/To Workspace2' */

  struct {
    void *LoggedData;
  } ToWorkspace20_PWORK;               /* '<S12>/To Workspace20' */

  struct {
    void *LoggedData;
  } ToWorkspace21_PWORK;               /* '<S12>/To Workspace21' */

  struct {
    void *LoggedData;
  } ToWorkspace22_PWORK;               /* '<S12>/To Workspace22' */

  struct {
    void *LoggedData;
  } ToWorkspace23_PWORK;               /* '<S12>/To Workspace23' */

  struct {
    void *LoggedData;
  } ToWorkspace24_PWORK;               /* '<S12>/To Workspace24' */

  struct {
    void *LoggedData;
  } ToWorkspace25_PWORK;               /* '<S12>/To Workspace25' */

  struct {
    void *LoggedData;
  } ToWorkspace26_PWORK;               /* '<S12>/To Workspace26' */

  struct {
    void *LoggedData;
  } ToWorkspace27_PWORK;               /* '<S12>/To Workspace27' */

  struct {
    void *LoggedData;
  } ToWorkspace28_PWORK;               /* '<S12>/To Workspace28' */

  struct {
    void *LoggedData;
  } ToWorkspace29_PWORK;               /* '<S12>/To Workspace29' */

  struct {
    void *LoggedData;
  } ToWorkspace3_PWORK;                /* '<S12>/To Workspace3' */

  struct {
    void *LoggedData;
  } ToWorkspace4_PWORK;                /* '<S12>/To Workspace4' */

  struct {
    void *LoggedData;
  } ToWorkspace5_PWORK;                /* '<S12>/To Workspace5' */

  struct {
    void *LoggedData;
  } ToWorkspace6_PWORK;                /* '<S12>/To Workspace6' */

  struct {
    void *LoggedData;
  } ToWorkspace7_PWORK;                /* '<S12>/To Workspace7' */

  struct {
    void *LoggedData;
  } ToWorkspace8_PWORK;                /* '<S12>/To Workspace8' */

  struct {
    void *LoggedData;
  } ToWorkspace9_PWORK;                /* '<S12>/To Workspace9' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_g;               /* '<S12>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S12>/To Workspace1' */

  struct {
    void *LoggedData;
  } ToWorkspace30_PWORK;               /* '<S12>/To Workspace30' */

  struct {
    void *LoggedData;
  } ToWorkspace32_PWORK;               /* '<S12>/To Workspace32' */

  struct {
    void *LoggedData;
  } ToWorkspace33_PWORK;               /* '<S12>/To Workspace33' */

  struct {
    void *LoggedData;
  } ToWorkspace31_PWORK;               /* '<S12>/To Workspace31' */

  robotics_slcore_internal_bloc_T obj_c;
                             /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T Memory_PreviousInput;      /* '<S145>/Memory' */
  boolean_T Memory_PreviousInput_d;    /* '<S91>/Memory' */
  boolean_T objisempty;                /* '<S11>/SourceBlock' */
  boolean_T objisempty_a;              /* '<S10>/SourceBlock' */
  boolean_T objisempty_n;              /* '<S9>/SourceBlock' */
  boolean_T objisempty_d;    /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_f;              /* '<S2>/ServiceCaller' */
};

/* Continuous states (default storage) */
struct X_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_p;          /* '<S101>/Integrator' */
  real_T Filter_CSTATE;                /* '<S96>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S207>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S202>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S47>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S42>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S155>/Integrator' */
  real_T Filter_CSTATE_fy;             /* '<S150>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S12>/Integrator1' */
};

/* State derivatives (default storage) */
struct XDot_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_p;          /* '<S101>/Integrator' */
  real_T Filter_CSTATE;                /* '<S96>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S207>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S202>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S47>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S42>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S155>/Integrator' */
  real_T Filter_CSTATE_fy;             /* '<S150>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S12>/Integrator1' */
};

/* State disabled  */
struct XDis_pid_control_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S12>/Integrator' */
  boolean_T Integrator_CSTATE_p;       /* '<S101>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S96>/Filter' */
  boolean_T Integrator_CSTATE_d;       /* '<S207>/Integrator' */
  boolean_T Filter_CSTATE_f;           /* '<S202>/Filter' */
  boolean_T Integrator_CSTATE_m;       /* '<S47>/Integrator' */
  boolean_T Filter_CSTATE_g;           /* '<S42>/Filter' */
  boolean_T Integrator_CSTATE_b;       /* '<S155>/Integrator' */
  boolean_T Filter_CSTATE_fy;          /* '<S150>/Filter' */
  boolean_T Integrator1_CSTATE;        /* '<S12>/Integrator1' */
};

/* Invariant block signals (default storage) */
struct ConstB_pid_control_V1_T {
  real_T motorspeed;                   /* '<S12>/Gain2' */
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
struct ConstP_pid_control_V1_T {
  /* Expression: x0
   * Referenced by: '<S12>/Integrator'
   */
  real_T Integrator_IC[12];
};

/* Real-time Model Data Structure */
struct tag_RTM_pid_control_V1_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_pid_control_V1_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  XDis_pid_control_V1_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[21];
  real_T odeF[3][21];
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

extern const ConstB_pid_control_V1_T pid_control_V1_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_pid_control_V1_T pid_control_V1_ConstP;

/* Class declaration for model pid_control_V1 */
class pid_control_V1
{
  /* public data and function members */
 public:
  /* Real-Time Model get method */
  RT_MODEL_pid_control_V1_T * getRTM();

  /* model start function */
  void start();

  /* Initial conditions function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  void terminate();

  /* Constructor */
  pid_control_V1();

  /* Destructor */
  ~pid_control_V1();

  /* private data and function members */
 private:
  /* Block signals */
  B_pid_control_V1_T pid_control_V1_B;

  /* Block states */
  DW_pid_control_V1_T pid_control_V1_DW;

  /* Block continuous states */
  X_pid_control_V1_T pid_control_V1_X;

  /* Block Continuous state disabled vector */
  XDis_pid_control_V1_T pid_control_V1_XDis;

  /* private member function(s) for subsystem '<S9>/Enabled Subsystem'*/
  static void pid_contr_EnabledSubsystem_Init(B_EnabledSubsystem_pid_contro_T
    *localB);
  static void pid_control_V1_EnabledSubsystem(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Float64 *rtu_In1, B_EnabledSubsystem_pid_contro_T *localB);

  /* private member function(s) for subsystem '<Root>'*/
  void pid_cont_Subscriber_setupImpl_o(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_con_Subscriber_setupImpl_on(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_contro_Subscriber_setupImpl(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_con_ServiceCaller_setupImpl(const ros_slros2_internal_block_Ser_T
    *obj);

  /* Global mass matrix */

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void pid_control_V1_derivatives();

  /* Real-Time Model */
  RT_MODEL_pid_control_V1_T pid_control_V1_M;
};

extern volatile boolean_T stopRequested;
extern volatile boolean_T runModel;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/Display' : Unused code path elimination
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
 * '<Root>' : 'pid_control_V1'
 * '<S1>'   : 'pid_control_V1/Blank Message'
 * '<S2>'   : 'pid_control_V1/Call Service'
 * '<S3>'   : 'pid_control_V1/MATLAB Function'
 * '<S4>'   : 'pid_control_V1/MATLAB Function1'
 * '<S5>'   : 'pid_control_V1/PID ALERON'
 * '<S6>'   : 'pid_control_V1/PID ELEVATOR'
 * '<S7>'   : 'pid_control_V1/PID MOTORS ALTURA'
 * '<S8>'   : 'pid_control_V1/PID TIIMON'
 * '<S9>'   : 'pid_control_V1/Subscribe'
 * '<S10>'  : 'pid_control_V1/Subscribe1'
 * '<S11>'  : 'pid_control_V1/Subscribe2'
 * '<S12>'  : 'pid_control_V1/Subsystem'
 * '<S13>'  : 'pid_control_V1/PID ALERON/Anti-windup'
 * '<S14>'  : 'pid_control_V1/PID ALERON/D Gain'
 * '<S15>'  : 'pid_control_V1/PID ALERON/External Derivative'
 * '<S16>'  : 'pid_control_V1/PID ALERON/Filter'
 * '<S17>'  : 'pid_control_V1/PID ALERON/Filter ICs'
 * '<S18>'  : 'pid_control_V1/PID ALERON/I Gain'
 * '<S19>'  : 'pid_control_V1/PID ALERON/Ideal P Gain'
 * '<S20>'  : 'pid_control_V1/PID ALERON/Ideal P Gain Fdbk'
 * '<S21>'  : 'pid_control_V1/PID ALERON/Integrator'
 * '<S22>'  : 'pid_control_V1/PID ALERON/Integrator ICs'
 * '<S23>'  : 'pid_control_V1/PID ALERON/N Copy'
 * '<S24>'  : 'pid_control_V1/PID ALERON/N Gain'
 * '<S25>'  : 'pid_control_V1/PID ALERON/P Copy'
 * '<S26>'  : 'pid_control_V1/PID ALERON/Parallel P Gain'
 * '<S27>'  : 'pid_control_V1/PID ALERON/Reset Signal'
 * '<S28>'  : 'pid_control_V1/PID ALERON/Saturation'
 * '<S29>'  : 'pid_control_V1/PID ALERON/Saturation Fdbk'
 * '<S30>'  : 'pid_control_V1/PID ALERON/Sum'
 * '<S31>'  : 'pid_control_V1/PID ALERON/Sum Fdbk'
 * '<S32>'  : 'pid_control_V1/PID ALERON/Tracking Mode'
 * '<S33>'  : 'pid_control_V1/PID ALERON/Tracking Mode Sum'
 * '<S34>'  : 'pid_control_V1/PID ALERON/Tsamp - Integral'
 * '<S35>'  : 'pid_control_V1/PID ALERON/Tsamp - Ngain'
 * '<S36>'  : 'pid_control_V1/PID ALERON/postSat Signal'
 * '<S37>'  : 'pid_control_V1/PID ALERON/preInt Signal'
 * '<S38>'  : 'pid_control_V1/PID ALERON/preSat Signal'
 * '<S39>'  : 'pid_control_V1/PID ALERON/Anti-windup/Passthrough'
 * '<S40>'  : 'pid_control_V1/PID ALERON/D Gain/Internal Parameters'
 * '<S41>'  : 'pid_control_V1/PID ALERON/External Derivative/Error'
 * '<S42>'  : 'pid_control_V1/PID ALERON/Filter/Cont. Filter'
 * '<S43>'  : 'pid_control_V1/PID ALERON/Filter ICs/Internal IC - Filter'
 * '<S44>'  : 'pid_control_V1/PID ALERON/I Gain/Internal Parameters'
 * '<S45>'  : 'pid_control_V1/PID ALERON/Ideal P Gain/Passthrough'
 * '<S46>'  : 'pid_control_V1/PID ALERON/Ideal P Gain Fdbk/Disabled'
 * '<S47>'  : 'pid_control_V1/PID ALERON/Integrator/Continuous'
 * '<S48>'  : 'pid_control_V1/PID ALERON/Integrator ICs/Internal IC'
 * '<S49>'  : 'pid_control_V1/PID ALERON/N Copy/Disabled'
 * '<S50>'  : 'pid_control_V1/PID ALERON/N Gain/Internal Parameters'
 * '<S51>'  : 'pid_control_V1/PID ALERON/P Copy/Disabled'
 * '<S52>'  : 'pid_control_V1/PID ALERON/Parallel P Gain/Internal Parameters'
 * '<S53>'  : 'pid_control_V1/PID ALERON/Reset Signal/Disabled'
 * '<S54>'  : 'pid_control_V1/PID ALERON/Saturation/Enabled'
 * '<S55>'  : 'pid_control_V1/PID ALERON/Saturation Fdbk/Disabled'
 * '<S56>'  : 'pid_control_V1/PID ALERON/Sum/Sum_PID'
 * '<S57>'  : 'pid_control_V1/PID ALERON/Sum Fdbk/Disabled'
 * '<S58>'  : 'pid_control_V1/PID ALERON/Tracking Mode/Disabled'
 * '<S59>'  : 'pid_control_V1/PID ALERON/Tracking Mode Sum/Passthrough'
 * '<S60>'  : 'pid_control_V1/PID ALERON/Tsamp - Integral/TsSignalSpecification'
 * '<S61>'  : 'pid_control_V1/PID ALERON/Tsamp - Ngain/Passthrough'
 * '<S62>'  : 'pid_control_V1/PID ALERON/postSat Signal/Forward_Path'
 * '<S63>'  : 'pid_control_V1/PID ALERON/preInt Signal/Internal PreInt'
 * '<S64>'  : 'pid_control_V1/PID ALERON/preSat Signal/Forward_Path'
 * '<S65>'  : 'pid_control_V1/PID ELEVATOR/Anti-windup'
 * '<S66>'  : 'pid_control_V1/PID ELEVATOR/D Gain'
 * '<S67>'  : 'pid_control_V1/PID ELEVATOR/External Derivative'
 * '<S68>'  : 'pid_control_V1/PID ELEVATOR/Filter'
 * '<S69>'  : 'pid_control_V1/PID ELEVATOR/Filter ICs'
 * '<S70>'  : 'pid_control_V1/PID ELEVATOR/I Gain'
 * '<S71>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain'
 * '<S72>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain Fdbk'
 * '<S73>'  : 'pid_control_V1/PID ELEVATOR/Integrator'
 * '<S74>'  : 'pid_control_V1/PID ELEVATOR/Integrator ICs'
 * '<S75>'  : 'pid_control_V1/PID ELEVATOR/N Copy'
 * '<S76>'  : 'pid_control_V1/PID ELEVATOR/N Gain'
 * '<S77>'  : 'pid_control_V1/PID ELEVATOR/P Copy'
 * '<S78>'  : 'pid_control_V1/PID ELEVATOR/Parallel P Gain'
 * '<S79>'  : 'pid_control_V1/PID ELEVATOR/Reset Signal'
 * '<S80>'  : 'pid_control_V1/PID ELEVATOR/Saturation'
 * '<S81>'  : 'pid_control_V1/PID ELEVATOR/Saturation Fdbk'
 * '<S82>'  : 'pid_control_V1/PID ELEVATOR/Sum'
 * '<S83>'  : 'pid_control_V1/PID ELEVATOR/Sum Fdbk'
 * '<S84>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode'
 * '<S85>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode Sum'
 * '<S86>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Integral'
 * '<S87>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Ngain'
 * '<S88>'  : 'pid_control_V1/PID ELEVATOR/postSat Signal'
 * '<S89>'  : 'pid_control_V1/PID ELEVATOR/preInt Signal'
 * '<S90>'  : 'pid_control_V1/PID ELEVATOR/preSat Signal'
 * '<S91>'  : 'pid_control_V1/PID ELEVATOR/Anti-windup/Cont. Clamping Parallel'
 * '<S92>'  : 'pid_control_V1/PID ELEVATOR/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S93>'  : 'pid_control_V1/PID ELEVATOR/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S94>'  : 'pid_control_V1/PID ELEVATOR/D Gain/Internal Parameters'
 * '<S95>'  : 'pid_control_V1/PID ELEVATOR/External Derivative/Error'
 * '<S96>'  : 'pid_control_V1/PID ELEVATOR/Filter/Cont. Filter'
 * '<S97>'  : 'pid_control_V1/PID ELEVATOR/Filter ICs/Internal IC - Filter'
 * '<S98>'  : 'pid_control_V1/PID ELEVATOR/I Gain/Internal Parameters'
 * '<S99>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain/Passthrough'
 * '<S100>' : 'pid_control_V1/PID ELEVATOR/Ideal P Gain Fdbk/Disabled'
 * '<S101>' : 'pid_control_V1/PID ELEVATOR/Integrator/Continuous'
 * '<S102>' : 'pid_control_V1/PID ELEVATOR/Integrator ICs/Internal IC'
 * '<S103>' : 'pid_control_V1/PID ELEVATOR/N Copy/Disabled'
 * '<S104>' : 'pid_control_V1/PID ELEVATOR/N Gain/Internal Parameters'
 * '<S105>' : 'pid_control_V1/PID ELEVATOR/P Copy/Disabled'
 * '<S106>' : 'pid_control_V1/PID ELEVATOR/Parallel P Gain/Internal Parameters'
 * '<S107>' : 'pid_control_V1/PID ELEVATOR/Reset Signal/Disabled'
 * '<S108>' : 'pid_control_V1/PID ELEVATOR/Saturation/Enabled'
 * '<S109>' : 'pid_control_V1/PID ELEVATOR/Saturation Fdbk/Disabled'
 * '<S110>' : 'pid_control_V1/PID ELEVATOR/Sum/Sum_PID'
 * '<S111>' : 'pid_control_V1/PID ELEVATOR/Sum Fdbk/Disabled'
 * '<S112>' : 'pid_control_V1/PID ELEVATOR/Tracking Mode/Disabled'
 * '<S113>' : 'pid_control_V1/PID ELEVATOR/Tracking Mode Sum/Passthrough'
 * '<S114>' : 'pid_control_V1/PID ELEVATOR/Tsamp - Integral/TsSignalSpecification'
 * '<S115>' : 'pid_control_V1/PID ELEVATOR/Tsamp - Ngain/Passthrough'
 * '<S116>' : 'pid_control_V1/PID ELEVATOR/postSat Signal/Forward_Path'
 * '<S117>' : 'pid_control_V1/PID ELEVATOR/preInt Signal/Internal PreInt'
 * '<S118>' : 'pid_control_V1/PID ELEVATOR/preSat Signal/Forward_Path'
 * '<S119>' : 'pid_control_V1/PID MOTORS ALTURA/Anti-windup'
 * '<S120>' : 'pid_control_V1/PID MOTORS ALTURA/D Gain'
 * '<S121>' : 'pid_control_V1/PID MOTORS ALTURA/External Derivative'
 * '<S122>' : 'pid_control_V1/PID MOTORS ALTURA/Filter'
 * '<S123>' : 'pid_control_V1/PID MOTORS ALTURA/Filter ICs'
 * '<S124>' : 'pid_control_V1/PID MOTORS ALTURA/I Gain'
 * '<S125>' : 'pid_control_V1/PID MOTORS ALTURA/Ideal P Gain'
 * '<S126>' : 'pid_control_V1/PID MOTORS ALTURA/Ideal P Gain Fdbk'
 * '<S127>' : 'pid_control_V1/PID MOTORS ALTURA/Integrator'
 * '<S128>' : 'pid_control_V1/PID MOTORS ALTURA/Integrator ICs'
 * '<S129>' : 'pid_control_V1/PID MOTORS ALTURA/N Copy'
 * '<S130>' : 'pid_control_V1/PID MOTORS ALTURA/N Gain'
 * '<S131>' : 'pid_control_V1/PID MOTORS ALTURA/P Copy'
 * '<S132>' : 'pid_control_V1/PID MOTORS ALTURA/Parallel P Gain'
 * '<S133>' : 'pid_control_V1/PID MOTORS ALTURA/Reset Signal'
 * '<S134>' : 'pid_control_V1/PID MOTORS ALTURA/Saturation'
 * '<S135>' : 'pid_control_V1/PID MOTORS ALTURA/Saturation Fdbk'
 * '<S136>' : 'pid_control_V1/PID MOTORS ALTURA/Sum'
 * '<S137>' : 'pid_control_V1/PID MOTORS ALTURA/Sum Fdbk'
 * '<S138>' : 'pid_control_V1/PID MOTORS ALTURA/Tracking Mode'
 * '<S139>' : 'pid_control_V1/PID MOTORS ALTURA/Tracking Mode Sum'
 * '<S140>' : 'pid_control_V1/PID MOTORS ALTURA/Tsamp - Integral'
 * '<S141>' : 'pid_control_V1/PID MOTORS ALTURA/Tsamp - Ngain'
 * '<S142>' : 'pid_control_V1/PID MOTORS ALTURA/postSat Signal'
 * '<S143>' : 'pid_control_V1/PID MOTORS ALTURA/preInt Signal'
 * '<S144>' : 'pid_control_V1/PID MOTORS ALTURA/preSat Signal'
 * '<S145>' : 'pid_control_V1/PID MOTORS ALTURA/Anti-windup/Cont. Clamping Parallel'
 * '<S146>' : 'pid_control_V1/PID MOTORS ALTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S147>' : 'pid_control_V1/PID MOTORS ALTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S148>' : 'pid_control_V1/PID MOTORS ALTURA/D Gain/Internal Parameters'
 * '<S149>' : 'pid_control_V1/PID MOTORS ALTURA/External Derivative/Error'
 * '<S150>' : 'pid_control_V1/PID MOTORS ALTURA/Filter/Cont. Filter'
 * '<S151>' : 'pid_control_V1/PID MOTORS ALTURA/Filter ICs/Internal IC - Filter'
 * '<S152>' : 'pid_control_V1/PID MOTORS ALTURA/I Gain/Internal Parameters'
 * '<S153>' : 'pid_control_V1/PID MOTORS ALTURA/Ideal P Gain/Passthrough'
 * '<S154>' : 'pid_control_V1/PID MOTORS ALTURA/Ideal P Gain Fdbk/Disabled'
 * '<S155>' : 'pid_control_V1/PID MOTORS ALTURA/Integrator/Continuous'
 * '<S156>' : 'pid_control_V1/PID MOTORS ALTURA/Integrator ICs/Internal IC'
 * '<S157>' : 'pid_control_V1/PID MOTORS ALTURA/N Copy/Disabled'
 * '<S158>' : 'pid_control_V1/PID MOTORS ALTURA/N Gain/Internal Parameters'
 * '<S159>' : 'pid_control_V1/PID MOTORS ALTURA/P Copy/Disabled'
 * '<S160>' : 'pid_control_V1/PID MOTORS ALTURA/Parallel P Gain/Internal Parameters'
 * '<S161>' : 'pid_control_V1/PID MOTORS ALTURA/Reset Signal/Disabled'
 * '<S162>' : 'pid_control_V1/PID MOTORS ALTURA/Saturation/Enabled'
 * '<S163>' : 'pid_control_V1/PID MOTORS ALTURA/Saturation Fdbk/Disabled'
 * '<S164>' : 'pid_control_V1/PID MOTORS ALTURA/Sum/Sum_PID'
 * '<S165>' : 'pid_control_V1/PID MOTORS ALTURA/Sum Fdbk/Disabled'
 * '<S166>' : 'pid_control_V1/PID MOTORS ALTURA/Tracking Mode/Disabled'
 * '<S167>' : 'pid_control_V1/PID MOTORS ALTURA/Tracking Mode Sum/Passthrough'
 * '<S168>' : 'pid_control_V1/PID MOTORS ALTURA/Tsamp - Integral/TsSignalSpecification'
 * '<S169>' : 'pid_control_V1/PID MOTORS ALTURA/Tsamp - Ngain/Passthrough'
 * '<S170>' : 'pid_control_V1/PID MOTORS ALTURA/postSat Signal/Forward_Path'
 * '<S171>' : 'pid_control_V1/PID MOTORS ALTURA/preInt Signal/Internal PreInt'
 * '<S172>' : 'pid_control_V1/PID MOTORS ALTURA/preSat Signal/Forward_Path'
 * '<S173>' : 'pid_control_V1/PID TIIMON/Anti-windup'
 * '<S174>' : 'pid_control_V1/PID TIIMON/D Gain'
 * '<S175>' : 'pid_control_V1/PID TIIMON/External Derivative'
 * '<S176>' : 'pid_control_V1/PID TIIMON/Filter'
 * '<S177>' : 'pid_control_V1/PID TIIMON/Filter ICs'
 * '<S178>' : 'pid_control_V1/PID TIIMON/I Gain'
 * '<S179>' : 'pid_control_V1/PID TIIMON/Ideal P Gain'
 * '<S180>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk'
 * '<S181>' : 'pid_control_V1/PID TIIMON/Integrator'
 * '<S182>' : 'pid_control_V1/PID TIIMON/Integrator ICs'
 * '<S183>' : 'pid_control_V1/PID TIIMON/N Copy'
 * '<S184>' : 'pid_control_V1/PID TIIMON/N Gain'
 * '<S185>' : 'pid_control_V1/PID TIIMON/P Copy'
 * '<S186>' : 'pid_control_V1/PID TIIMON/Parallel P Gain'
 * '<S187>' : 'pid_control_V1/PID TIIMON/Reset Signal'
 * '<S188>' : 'pid_control_V1/PID TIIMON/Saturation'
 * '<S189>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk'
 * '<S190>' : 'pid_control_V1/PID TIIMON/Sum'
 * '<S191>' : 'pid_control_V1/PID TIIMON/Sum Fdbk'
 * '<S192>' : 'pid_control_V1/PID TIIMON/Tracking Mode'
 * '<S193>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum'
 * '<S194>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral'
 * '<S195>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain'
 * '<S196>' : 'pid_control_V1/PID TIIMON/postSat Signal'
 * '<S197>' : 'pid_control_V1/PID TIIMON/preInt Signal'
 * '<S198>' : 'pid_control_V1/PID TIIMON/preSat Signal'
 * '<S199>' : 'pid_control_V1/PID TIIMON/Anti-windup/Passthrough'
 * '<S200>' : 'pid_control_V1/PID TIIMON/D Gain/Internal Parameters'
 * '<S201>' : 'pid_control_V1/PID TIIMON/External Derivative/Error'
 * '<S202>' : 'pid_control_V1/PID TIIMON/Filter/Cont. Filter'
 * '<S203>' : 'pid_control_V1/PID TIIMON/Filter ICs/Internal IC - Filter'
 * '<S204>' : 'pid_control_V1/PID TIIMON/I Gain/Internal Parameters'
 * '<S205>' : 'pid_control_V1/PID TIIMON/Ideal P Gain/Passthrough'
 * '<S206>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk/Disabled'
 * '<S207>' : 'pid_control_V1/PID TIIMON/Integrator/Continuous'
 * '<S208>' : 'pid_control_V1/PID TIIMON/Integrator ICs/Internal IC'
 * '<S209>' : 'pid_control_V1/PID TIIMON/N Copy/Disabled'
 * '<S210>' : 'pid_control_V1/PID TIIMON/N Gain/Internal Parameters'
 * '<S211>' : 'pid_control_V1/PID TIIMON/P Copy/Disabled'
 * '<S212>' : 'pid_control_V1/PID TIIMON/Parallel P Gain/Internal Parameters'
 * '<S213>' : 'pid_control_V1/PID TIIMON/Reset Signal/Disabled'
 * '<S214>' : 'pid_control_V1/PID TIIMON/Saturation/Enabled'
 * '<S215>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk/Disabled'
 * '<S216>' : 'pid_control_V1/PID TIIMON/Sum/Sum_PID'
 * '<S217>' : 'pid_control_V1/PID TIIMON/Sum Fdbk/Disabled'
 * '<S218>' : 'pid_control_V1/PID TIIMON/Tracking Mode/Disabled'
 * '<S219>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum/Passthrough'
 * '<S220>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral/TsSignalSpecification'
 * '<S221>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain/Passthrough'
 * '<S222>' : 'pid_control_V1/PID TIIMON/postSat Signal/Forward_Path'
 * '<S223>' : 'pid_control_V1/PID TIIMON/preInt Signal/Internal PreInt'
 * '<S224>' : 'pid_control_V1/PID TIIMON/preSat Signal/Forward_Path'
 * '<S225>' : 'pid_control_V1/Subscribe/Enabled Subsystem'
 * '<S226>' : 'pid_control_V1/Subscribe1/Enabled Subsystem'
 * '<S227>' : 'pid_control_V1/Subscribe2/Enabled Subsystem'
 * '<S228>' : 'pid_control_V1/Subsystem/MATLAB Function'
 */
#endif                                 /* pid_control_V1_h_ */
