/*
 * pid_control_V1.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.30
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Sat Feb 28 18:41:08 2026
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

/* Block signals (default storage) */
struct B_pid_control_V1_T {
  SL_Bus_gazebo_msgs_SetEntityStateRequest BusAssignment;/* '<Root>/Bus Assignment' */
  real_T x[12];                        /* '<S7>/Integrator' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function1' */
  uint8_T stringOut_l[128];            /* '<Root>/MATLAB Function' */
  real_T FA_b_tmp[9];
  real_T control_vector[5];
  char_T b_zeroDelimTopic[25];
  real_T wbe_b[3];
  real_T FE2_b[3];
  real_T F_b[3];
  real_T Q[3];
  real_T dv[2];
  real_T FilterCoefficient;            /* '<S45>/Filter Coefficient' */
  real_T Sum;                          /* '<S51>/Sum' */
  real_T FilterCoefficient_g;          /* '<S99>/Filter Coefficient' */
  real_T Saturation;                   /* '<S103>/Saturation' */
  real_T Gain5;                        /* '<Root>/Gain5' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T Power;                        /* '<S7>/Product2' */
  real_T Gain3;                        /* '<S7>/Gain3' */
  real_T powerdemand;                  /* '<S7>/Divide' */
  real_T loadtorque;                   /* '<S7>/Divide1' */
  real_T Switch;                       /* '<S86>/Switch' */
  real_T IntegralGain;                 /* '<S39>/Integral Gain' */
  real_T EnergykWh;                    /* '<S7>/Gain1' */
  real_T XDOT[40];                     /* '<S7>/MATLAB Function' */
  real_T u2;
  real_T Va;
  real_T alpha;
  real_T beta;
  real_T hw;
  real_T Q_m;
  real_T CL_w_OGE;
  real_T CL_h_OGE;
  real_T CL_h_IGE;
  real_T CD_iw_IGE;
  real_T CD_ih_IGE;
  real_T CQ;
  real_T Cl;
  real_T Cn;
  real_T Vd1;
  real_T Vd2;
  real_T SignPreSat;                   /* '<S86>/SignPreSat' */
  real_T ZeroGain;                     /* '<S86>/ZeroGain' */
  real_T Fg_b_tmp;
  real_T Fg_b_tmp_c;
  real_T Fg_b_tmp_k;
  real_T FE2_b_tmp;
  real_T Fg_b_tmp_cx;
  real_T FE2_b_tmp_b;
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
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function1' */
  uint32_T lengthOut_e;                /* '<Root>/MATLAB Function' */
  boolean_T AND3;                      /* '<S86>/AND3' */
  boolean_T Memory;                    /* '<S86>/Memory' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_pid_control_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<Root>/To Workspace' */

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
  } ToWorkspace_PWORK_g;               /* '<S7>/To Workspace' */

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
  boolean_T Memory_PreviousInput;      /* '<S86>/Memory' */
  boolean_T objisempty;      /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_f;              /* '<S2>/ServiceCaller' */
};

/* Continuous states (default storage) */
struct X_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S7>/Integrator' */
  real_T Integrator_CSTATE_p;          /* '<S42>/Integrator' */
  real_T Filter_CSTATE;                /* '<S37>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S96>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S91>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S7>/Integrator1' */
};

/* State derivatives (default storage) */
struct XDot_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S7>/Integrator' */
  real_T Integrator_CSTATE_p;          /* '<S42>/Integrator' */
  real_T Filter_CSTATE;                /* '<S37>/Filter' */
  real_T Integrator_CSTATE_b;          /* '<S96>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S91>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S7>/Integrator1' */
};

/* State disabled  */
struct XDis_pid_control_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S7>/Integrator' */
  boolean_T Integrator_CSTATE_p;       /* '<S42>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S37>/Filter' */
  boolean_T Integrator_CSTATE_b;       /* '<S96>/Integrator' */
  boolean_T Filter_CSTATE_f;           /* '<S91>/Filter' */
  boolean_T Integrator1_CSTATE;        /* '<S7>/Integrator1' */
};

/* Invariant block signals (default storage) */
struct ConstB_pid_control_V1_T {
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
struct ConstP_pid_control_V1_T {
  /* Expression: x0
   * Referenced by: '<S7>/Integrator'
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

  /* private member function(s) for subsystem '<Root>'*/
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
 * '<S5>'   : 'pid_control_V1/PID ELEVATOR'
 * '<S6>'   : 'pid_control_V1/PID MOTORS AlTURA'
 * '<S7>'   : 'pid_control_V1/Subsystem'
 * '<S8>'   : 'pid_control_V1/PID ELEVATOR/Anti-windup'
 * '<S9>'   : 'pid_control_V1/PID ELEVATOR/D Gain'
 * '<S10>'  : 'pid_control_V1/PID ELEVATOR/External Derivative'
 * '<S11>'  : 'pid_control_V1/PID ELEVATOR/Filter'
 * '<S12>'  : 'pid_control_V1/PID ELEVATOR/Filter ICs'
 * '<S13>'  : 'pid_control_V1/PID ELEVATOR/I Gain'
 * '<S14>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain'
 * '<S15>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain Fdbk'
 * '<S16>'  : 'pid_control_V1/PID ELEVATOR/Integrator'
 * '<S17>'  : 'pid_control_V1/PID ELEVATOR/Integrator ICs'
 * '<S18>'  : 'pid_control_V1/PID ELEVATOR/N Copy'
 * '<S19>'  : 'pid_control_V1/PID ELEVATOR/N Gain'
 * '<S20>'  : 'pid_control_V1/PID ELEVATOR/P Copy'
 * '<S21>'  : 'pid_control_V1/PID ELEVATOR/Parallel P Gain'
 * '<S22>'  : 'pid_control_V1/PID ELEVATOR/Reset Signal'
 * '<S23>'  : 'pid_control_V1/PID ELEVATOR/Saturation'
 * '<S24>'  : 'pid_control_V1/PID ELEVATOR/Saturation Fdbk'
 * '<S25>'  : 'pid_control_V1/PID ELEVATOR/Sum'
 * '<S26>'  : 'pid_control_V1/PID ELEVATOR/Sum Fdbk'
 * '<S27>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode'
 * '<S28>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode Sum'
 * '<S29>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Integral'
 * '<S30>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Ngain'
 * '<S31>'  : 'pid_control_V1/PID ELEVATOR/postSat Signal'
 * '<S32>'  : 'pid_control_V1/PID ELEVATOR/preInt Signal'
 * '<S33>'  : 'pid_control_V1/PID ELEVATOR/preSat Signal'
 * '<S34>'  : 'pid_control_V1/PID ELEVATOR/Anti-windup/Passthrough'
 * '<S35>'  : 'pid_control_V1/PID ELEVATOR/D Gain/Internal Parameters'
 * '<S36>'  : 'pid_control_V1/PID ELEVATOR/External Derivative/Error'
 * '<S37>'  : 'pid_control_V1/PID ELEVATOR/Filter/Cont. Filter'
 * '<S38>'  : 'pid_control_V1/PID ELEVATOR/Filter ICs/Internal IC - Filter'
 * '<S39>'  : 'pid_control_V1/PID ELEVATOR/I Gain/Internal Parameters'
 * '<S40>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain/Passthrough'
 * '<S41>'  : 'pid_control_V1/PID ELEVATOR/Ideal P Gain Fdbk/Disabled'
 * '<S42>'  : 'pid_control_V1/PID ELEVATOR/Integrator/Continuous'
 * '<S43>'  : 'pid_control_V1/PID ELEVATOR/Integrator ICs/Internal IC'
 * '<S44>'  : 'pid_control_V1/PID ELEVATOR/N Copy/Disabled'
 * '<S45>'  : 'pid_control_V1/PID ELEVATOR/N Gain/Internal Parameters'
 * '<S46>'  : 'pid_control_V1/PID ELEVATOR/P Copy/Disabled'
 * '<S47>'  : 'pid_control_V1/PID ELEVATOR/Parallel P Gain/Internal Parameters'
 * '<S48>'  : 'pid_control_V1/PID ELEVATOR/Reset Signal/Disabled'
 * '<S49>'  : 'pid_control_V1/PID ELEVATOR/Saturation/Passthrough'
 * '<S50>'  : 'pid_control_V1/PID ELEVATOR/Saturation Fdbk/Disabled'
 * '<S51>'  : 'pid_control_V1/PID ELEVATOR/Sum/Sum_PID'
 * '<S52>'  : 'pid_control_V1/PID ELEVATOR/Sum Fdbk/Disabled'
 * '<S53>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode/Disabled'
 * '<S54>'  : 'pid_control_V1/PID ELEVATOR/Tracking Mode Sum/Passthrough'
 * '<S55>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Integral/TsSignalSpecification'
 * '<S56>'  : 'pid_control_V1/PID ELEVATOR/Tsamp - Ngain/Passthrough'
 * '<S57>'  : 'pid_control_V1/PID ELEVATOR/postSat Signal/Forward_Path'
 * '<S58>'  : 'pid_control_V1/PID ELEVATOR/preInt Signal/Internal PreInt'
 * '<S59>'  : 'pid_control_V1/PID ELEVATOR/preSat Signal/Forward_Path'
 * '<S60>'  : 'pid_control_V1/PID MOTORS AlTURA/Anti-windup'
 * '<S61>'  : 'pid_control_V1/PID MOTORS AlTURA/D Gain'
 * '<S62>'  : 'pid_control_V1/PID MOTORS AlTURA/External Derivative'
 * '<S63>'  : 'pid_control_V1/PID MOTORS AlTURA/Filter'
 * '<S64>'  : 'pid_control_V1/PID MOTORS AlTURA/Filter ICs'
 * '<S65>'  : 'pid_control_V1/PID MOTORS AlTURA/I Gain'
 * '<S66>'  : 'pid_control_V1/PID MOTORS AlTURA/Ideal P Gain'
 * '<S67>'  : 'pid_control_V1/PID MOTORS AlTURA/Ideal P Gain Fdbk'
 * '<S68>'  : 'pid_control_V1/PID MOTORS AlTURA/Integrator'
 * '<S69>'  : 'pid_control_V1/PID MOTORS AlTURA/Integrator ICs'
 * '<S70>'  : 'pid_control_V1/PID MOTORS AlTURA/N Copy'
 * '<S71>'  : 'pid_control_V1/PID MOTORS AlTURA/N Gain'
 * '<S72>'  : 'pid_control_V1/PID MOTORS AlTURA/P Copy'
 * '<S73>'  : 'pid_control_V1/PID MOTORS AlTURA/Parallel P Gain'
 * '<S74>'  : 'pid_control_V1/PID MOTORS AlTURA/Reset Signal'
 * '<S75>'  : 'pid_control_V1/PID MOTORS AlTURA/Saturation'
 * '<S76>'  : 'pid_control_V1/PID MOTORS AlTURA/Saturation Fdbk'
 * '<S77>'  : 'pid_control_V1/PID MOTORS AlTURA/Sum'
 * '<S78>'  : 'pid_control_V1/PID MOTORS AlTURA/Sum Fdbk'
 * '<S79>'  : 'pid_control_V1/PID MOTORS AlTURA/Tracking Mode'
 * '<S80>'  : 'pid_control_V1/PID MOTORS AlTURA/Tracking Mode Sum'
 * '<S81>'  : 'pid_control_V1/PID MOTORS AlTURA/Tsamp - Integral'
 * '<S82>'  : 'pid_control_V1/PID MOTORS AlTURA/Tsamp - Ngain'
 * '<S83>'  : 'pid_control_V1/PID MOTORS AlTURA/postSat Signal'
 * '<S84>'  : 'pid_control_V1/PID MOTORS AlTURA/preInt Signal'
 * '<S85>'  : 'pid_control_V1/PID MOTORS AlTURA/preSat Signal'
 * '<S86>'  : 'pid_control_V1/PID MOTORS AlTURA/Anti-windup/Cont. Clamping Parallel'
 * '<S87>'  : 'pid_control_V1/PID MOTORS AlTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S88>'  : 'pid_control_V1/PID MOTORS AlTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S89>'  : 'pid_control_V1/PID MOTORS AlTURA/D Gain/Internal Parameters'
 * '<S90>'  : 'pid_control_V1/PID MOTORS AlTURA/External Derivative/Error'
 * '<S91>'  : 'pid_control_V1/PID MOTORS AlTURA/Filter/Cont. Filter'
 * '<S92>'  : 'pid_control_V1/PID MOTORS AlTURA/Filter ICs/Internal IC - Filter'
 * '<S93>'  : 'pid_control_V1/PID MOTORS AlTURA/I Gain/Internal Parameters'
 * '<S94>'  : 'pid_control_V1/PID MOTORS AlTURA/Ideal P Gain/Passthrough'
 * '<S95>'  : 'pid_control_V1/PID MOTORS AlTURA/Ideal P Gain Fdbk/Disabled'
 * '<S96>'  : 'pid_control_V1/PID MOTORS AlTURA/Integrator/Continuous'
 * '<S97>'  : 'pid_control_V1/PID MOTORS AlTURA/Integrator ICs/Internal IC'
 * '<S98>'  : 'pid_control_V1/PID MOTORS AlTURA/N Copy/Disabled'
 * '<S99>'  : 'pid_control_V1/PID MOTORS AlTURA/N Gain/Internal Parameters'
 * '<S100>' : 'pid_control_V1/PID MOTORS AlTURA/P Copy/Disabled'
 * '<S101>' : 'pid_control_V1/PID MOTORS AlTURA/Parallel P Gain/Internal Parameters'
 * '<S102>' : 'pid_control_V1/PID MOTORS AlTURA/Reset Signal/Disabled'
 * '<S103>' : 'pid_control_V1/PID MOTORS AlTURA/Saturation/Enabled'
 * '<S104>' : 'pid_control_V1/PID MOTORS AlTURA/Saturation Fdbk/Disabled'
 * '<S105>' : 'pid_control_V1/PID MOTORS AlTURA/Sum/Sum_PID'
 * '<S106>' : 'pid_control_V1/PID MOTORS AlTURA/Sum Fdbk/Disabled'
 * '<S107>' : 'pid_control_V1/PID MOTORS AlTURA/Tracking Mode/Disabled'
 * '<S108>' : 'pid_control_V1/PID MOTORS AlTURA/Tracking Mode Sum/Passthrough'
 * '<S109>' : 'pid_control_V1/PID MOTORS AlTURA/Tsamp - Integral/TsSignalSpecification'
 * '<S110>' : 'pid_control_V1/PID MOTORS AlTURA/Tsamp - Ngain/Passthrough'
 * '<S111>' : 'pid_control_V1/PID MOTORS AlTURA/postSat Signal/Forward_Path'
 * '<S112>' : 'pid_control_V1/PID MOTORS AlTURA/preInt Signal/Internal PreInt'
 * '<S113>' : 'pid_control_V1/PID MOTORS AlTURA/preSat Signal/Forward_Path'
 * '<S114>' : 'pid_control_V1/Subsystem/MATLAB Function'
 */
#endif                                 /* pid_control_V1_h_ */
