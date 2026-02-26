/*
 * open_loop_V1.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "open_loop_V1".
 *
 * Model version              : 12.16
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Thu Feb 26 12:25:33 2026
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
  real_T x[12];                        /* '<S4>/Integrator' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function' */
  real_T c_theta_g[9];
  char_T b_zeroDelimTopic[25];
  real_T Vb_w[3];
  real_T wbe_b[3];
  real_T F_b[3];
  real_T CD_iw_IGE_g[3];
  real_T Mcg_b[3];
  real_T dv[2];
  real_T Power;                        /* '<S4>/Product2' */
  real_T Gain3;                        /* '<S4>/Gain3' */
  real_T powerdemand;                  /* '<S4>/Divide' */
  real_T loadtorque;                   /* '<S4>/Divide1' */
  real_T EnergykWh;                    /* '<S4>/Gain1' */
  real_T XDOT[40];                     /* '<S4>/MATLAB Function' */
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
  real_T CL_h_IGE;
  real_T CD_iw_IGE;
  real_T CD_ih_IGE;
  real_T Cl;
  real_T Cm;
  real_T Cn;
  real_T Vd1;
  real_T FA_b_tmp;
  real_T FA_b_tmp_m;
  real_T Va_b_idx_0;
  real_T Va_b_idx_2;
  real_T Mcg_b_idx_1;
  real_T Mcg_b_idx_2;
  real_T FA_b_idx_0;
  real_T FA_b_idx_1;
  real_T FA_b_idx_2;
  real_T c_theta_tmp;
  real_T c_theta_tmp_c;
  real_T Vb_w_tmp;
  real_T c_theta_tmp_k;
  real_T c_theta_tmp_cx;
  real_T c_theta_tmp_b;
  real_T c_theta_tmp_p;
  real_T c_theta_tmp_cv;
  real_T d;
  real_T d1;
  int32_T i;
  int32_T i_f;
  int32_T i1;
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function' */
  boolean_T serverAvailableOnTime;
  boolean_T b;
  SL_Bus_gazebo_msgs_SetEntityStateResponse r;
};

/* Block states (default storage) for system '<Root>' */
struct DW_open_loop_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S4>/To Workspace1' */

  struct {
    void *LoggedData;
  } ToWorkspace10_PWORK;               /* '<S4>/To Workspace10' */

  struct {
    void *LoggedData;
  } ToWorkspace11_PWORK;               /* '<S4>/To Workspace11' */

  struct {
    void *LoggedData;
  } ToWorkspace12_PWORK;               /* '<S4>/To Workspace12' */

  struct {
    void *LoggedData;
  } ToWorkspace13_PWORK;               /* '<S4>/To Workspace13' */

  struct {
    void *LoggedData;
  } ToWorkspace14_PWORK;               /* '<S4>/To Workspace14' */

  struct {
    void *LoggedData;
  } ToWorkspace15_PWORK;               /* '<S4>/To Workspace15' */

  struct {
    void *LoggedData;
  } ToWorkspace16_PWORK;               /* '<S4>/To Workspace16' */

  struct {
    void *LoggedData;
  } ToWorkspace17_PWORK;               /* '<S4>/To Workspace17' */

  struct {
    void *LoggedData;
  } ToWorkspace18_PWORK;               /* '<S4>/To Workspace18' */

  struct {
    void *LoggedData;
  } ToWorkspace19_PWORK;               /* '<S4>/To Workspace19' */

  struct {
    void *LoggedData;
  } ToWorkspace2_PWORK;                /* '<S4>/To Workspace2' */

  struct {
    void *LoggedData;
  } ToWorkspace20_PWORK;               /* '<S4>/To Workspace20' */

  struct {
    void *LoggedData;
  } ToWorkspace21_PWORK;               /* '<S4>/To Workspace21' */

  struct {
    void *LoggedData;
  } ToWorkspace22_PWORK;               /* '<S4>/To Workspace22' */

  struct {
    void *LoggedData;
  } ToWorkspace23_PWORK;               /* '<S4>/To Workspace23' */

  struct {
    void *LoggedData;
  } ToWorkspace24_PWORK;               /* '<S4>/To Workspace24' */

  struct {
    void *LoggedData;
  } ToWorkspace25_PWORK;               /* '<S4>/To Workspace25' */

  struct {
    void *LoggedData;
  } ToWorkspace26_PWORK;               /* '<S4>/To Workspace26' */

  struct {
    void *LoggedData;
  } ToWorkspace27_PWORK;               /* '<S4>/To Workspace27' */

  struct {
    void *LoggedData;
  } ToWorkspace28_PWORK;               /* '<S4>/To Workspace28' */

  struct {
    void *LoggedData;
  } ToWorkspace29_PWORK;               /* '<S4>/To Workspace29' */

  struct {
    void *LoggedData;
  } ToWorkspace3_PWORK;                /* '<S4>/To Workspace3' */

  struct {
    void *LoggedData;
  } ToWorkspace4_PWORK;                /* '<S4>/To Workspace4' */

  struct {
    void *LoggedData;
  } ToWorkspace5_PWORK;                /* '<S4>/To Workspace5' */

  struct {
    void *LoggedData;
  } ToWorkspace6_PWORK;                /* '<S4>/To Workspace6' */

  struct {
    void *LoggedData;
  } ToWorkspace7_PWORK;                /* '<S4>/To Workspace7' */

  struct {
    void *LoggedData;
  } ToWorkspace8_PWORK;                /* '<S4>/To Workspace8' */

  struct {
    void *LoggedData;
  } ToWorkspace9_PWORK;                /* '<S4>/To Workspace9' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<S4>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace30_PWORK;               /* '<S4>/To Workspace30' */

  struct {
    void *LoggedData;
  } ToWorkspace32_PWORK;               /* '<S4>/To Workspace32' */

  struct {
    void *LoggedData;
  } ToWorkspace33_PWORK;               /* '<S4>/To Workspace33' */

  struct {
    void *LoggedData;
  } ToWorkspace31_PWORK;               /* '<S4>/To Workspace31' */

  robotics_slcore_internal_bloc_T obj_c;
                             /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty;      /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_f;              /* '<S2>/ServiceCaller' */
};

/* Continuous states (default storage) */
struct X_open_loop_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S4>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<S4>/Integrator1' */
};

/* State derivatives (default storage) */
struct XDot_open_loop_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S4>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<S4>/Integrator1' */
};

/* State disabled  */
struct XDis_open_loop_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S4>/Integrator' */
  boolean_T Integrator1_CSTATE;        /* '<S4>/Integrator1' */
};

/* Invariant block signals (default storage) */
struct ConstB_open_loop_V1_T {
  real_T motorspeed;                   /* '<S4>/Gain2' */
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
   * Referenced by: '<S4>/Integrator'
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
  real_T odeY[13];
  real_T odeF[3][13];
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
 * '<S4>'   : 'open_loop_V1/Subsystem'
 * '<S5>'   : 'open_loop_V1/Subsystem/MATLAB Function'
 */
#endif                                 /* open_loop_V1_h_ */
