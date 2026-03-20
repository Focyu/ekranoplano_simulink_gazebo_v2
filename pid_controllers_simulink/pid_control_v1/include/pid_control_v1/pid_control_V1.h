/*
 * pid_control_V1.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.93
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Fri Mar 20 12:02:13 2026
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
#include <string.h>

extern "C"
{

#include "rtGetInf.h"

}

extern "C"
{

#include "rtGetNaN.h"

}

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

/* Block signals for system '<S11>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_contro_T {
  SL_Bus_std_msgs_Float64 In1;         /* '<S283>/In1' */
};

/* Block signals for system '<S14>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_cont_p_T {
  SL_Bus_std_msgs_Int64 In1;           /* '<S286>/In1' */
};

/* Block signals for system '<S291>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_cont_n_T {
  SL_Bus_std_msgs_Bool In1;            /* '<S334>/In1' */
};

/* Block signals (default storage) */
struct B_pid_control_V1_T {
  SL_Bus_gazebo_msgs_SetEntityStateRequest BusAssignment;/* '<Root>/Bus Assignment' */
  real_T x[12];                        /* '<S16>/Integrator' */
  real_T RotationAnglestoDirectionCo[9];
                        /* '<S16>/Rotation Angles to Direction Cosine Matrix' */
  real_T FA_b_tmp[9];
  real_T TmpSignalConversionAtSFunct[5];/* '<S16>/MATLAB Function' */
  char_T b_zeroDelimTopic[25];
  real_T wbe_b[3];
  real_T FA_b[3];
  real_T F_b[3];
  real_T Product_be[3];                /* '<S319>/Product' */
  real_T wbe_b_m[3];
  char_T b_zeroDelimTopic_c[22];
  char_T b_zeroDelimTopic_k[22];
  char_T b_zeroDelimTopic_cx[21];
  char_T b_zeroDelimTopic_b[21];
  char_T b_zeroDelimTopic_p[17];
  char_T b_zeroDelimTopic_cv[17];
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF deadline_f;
  sJ4ih70VmKcvCeguWN0mNVF deadline_g;
  real_T w_n[2];                       /* '<S310>/w' */
  real_T w1_c[2];                      /* '<S310>/w1' */
  real_T w_d[2];                       /* '<S309>/w' */
  real_T w_e0[2];                      /* '<S308>/w' */
  real_T UnaryMinus[2];                /* '<S308>/Unary Minus' */
  real_T w_o[2];                       /* '<S307>/w' */
  real_T sigma_w[2];                   /* '<S307>/sigma_w' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function1' */
  uint8_T stringOut_l[128];            /* '<Root>/MATLAB Function' */
  real_T frac[2];
  real_T dv[2];
  real_T Switch3;                      /* '<Root>/Switch3' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T FilterCoefficient;            /* '<S110>/Filter Coefficient' */
  real_T Saturation;                   /* '<S114>/Saturation' */
  real_T Saturation1;                  /* '<Root>/Saturation1' */
  real_T FilterCoefficient_c;          /* '<S56>/Filter Coefficient' */
  real_T Saturation_k;                 /* '<S60>/Saturation' */
  real_T Saturation_i;                 /* '<Root>/Saturation' */
  real_T FilterCoefficient_m;          /* '<S162>/Filter Coefficient' */
  real_T Saturation_f;                 /* '<S166>/Saturation' */
  real_T Switch4;                      /* '<Root>/Switch4' */
  real_T Switch5;                      /* '<Root>/Switch5' */
  real_T Switch6;                      /* '<Root>/Switch6' */
  real_T Switch2;                      /* '<Root>/Switch2' */
  real_T Switch;                       /* '<Root>/Switch' */
  real_T FilterCoefficient_p;          /* '<S214>/Filter Coefficient' */
  real_T Saturation_m;                 /* '<S218>/Saturation' */
  real_T Switch_n;                     /* '<S43>/Switch' */
  real_T Switch_k;                     /* '<S97>/Switch' */
  real_T SumI4;                        /* '<S151>/SumI4' */
  real_T IntegralGain;                 /* '<S208>/Integral Gain' */
  real_T FilterCoefficient_cv;         /* '<S268>/Filter Coefficient' */
  real_T Switch_j;                     /* '<S255>/Switch' */
  real_T Saturation_o;                 /* '<S272>/Saturation' */
  real_T Memory[3];                    /* '<S16>/Memory' */
  real_T Memory1[3];                   /* '<S16>/Memory1' */
  real_T Power;                        /* '<S16>/Product2' */
  real_T Gain3;                        /* '<S16>/Gain3' */
  real_T EnergykWh;                    /* '<S16>/Gain1' */
  real_T powerdemand;                  /* '<S16>/Divide' */
  real_T loadtorque;                   /* '<S16>/Divide1' */
  real_T Output;                       /* '<S288>/Output' */
  real_T Product[4];                   /* '<S306>/Product' */
  real_T Sum[3];                       /* '<S16>/Sum' */
  real_T Sum1[3];                      /* '<S16>/Sum1' */
  real_T XDOT[40];                     /* '<S16>/MATLAB Function' */
  real_T w[2];                         /* '<S312>/w' */
  real_T w_a[2];                       /* '<S312>/w ' */
  real_T LwgV1[2];                     /* '<S312>/Lwg//V 1' */
  real_T w_g[2];                       /* '<S311>/w' */
  real_T w_e[2];                       /* '<S311>/w ' */
  real_T w1[2];                        /* '<S311>/w 1' */
  real_T chi;
  real_T u2;
  real_T Q;
  real_T Dtot;
  real_T Ltot;
  real_T CQ;
  real_T Cl;
  real_T Cn;
  real_T Vd1;
  real_T Tp1;
  real_T Tp2;
  real_T c_phi;
  real_T s_phi;
  real_T c_the;
  real_T s_the;
  real_T c_psi;
  real_T s_psi;
  real_T sina;
  real_T sinb;
  real_T sinc;
  real_T cosa;
  real_T cosb;
  real_T cosc;
  real_T SignPreSat;                   /* '<S97>/SignPreSat' */
  real_T SignPreSat_a;                 /* '<S43>/SignPreSat' */
  real_T Sum1_g;                       /* '<Root>/Sum1' */
  real_T Sum_hl;                       /* '<S168>/Sum' */
  real_T Sum5;                         /* '<Root>/Sum5' */
  real_T FE1_b_idx_1;
  real_T Mcg_b_idx_2;
  real_T Mcg_b_idx_0;
  real_T FE2_b_idx_0;
  real_T FE2_b_idx_2;
  real_T Fg_b_idx_2;
  real_T Fg_b_idx_1;
  real_T c_the_tmp;
  real_T c_the_tmp_g;
  SL_Bus_std_msgs_Float64 SourceBlock_o2_k;/* '<S294>/SourceBlock' */
  SL_Bus_std_msgs_Float64 SourceBlock_o2_p;/* '<S293>/SourceBlock' */
  SL_Bus_std_msgs_Float64 SourceBlock_o2_g2;/* '<S12>/SourceBlock' */
  SL_Bus_std_msgs_Int64 SourceBlock_o2_g;/* '<S14>/SourceBlock' */
  SL_Bus_std_msgs_Int64 SourceBlock_o2;/* '<S15>/SourceBlock' */
  uint32_T bpIndex[2];
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function1' */
  uint32_T lengthOut_e;                /* '<Root>/MATLAB Function' */
  boolean_T AND3;                      /* '<S43>/AND3' */
  boolean_T Memory_a;                  /* '<S43>/Memory' */
  boolean_T AND3_e;                    /* '<S97>/AND3' */
  boolean_T Memory_n;                  /* '<S97>/Memory' */
  boolean_T AND3_c;                    /* '<S255>/AND3' */
  boolean_T Memory_h;                  /* '<S255>/Memory' */
  boolean_T SourceBlock_o1;            /* '<S294>/SourceBlock' */
  boolean_T SourceBlock_o1_c;          /* '<S293>/SourceBlock' */
  boolean_T SourceBlock_o1_k;          /* '<S292>/SourceBlock' */
  boolean_T SourceBlock_o1_h;          /* '<S291>/SourceBlock' */
  boolean_T SourceBlock_o1_d;          /* '<S15>/SourceBlock' */
  boolean_T SourceBlock_o1_g;          /* '<S14>/SourceBlock' */
  boolean_T SourceBlock_o1_o;          /* '<S13>/SourceBlock' */
  boolean_T SourceBlock_o1_gx;         /* '<S12>/SourceBlock' */
  boolean_T SourceBlock_o1_l;          /* '<S11>/SourceBlock' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_pu;/* '<S294>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_k;/* '<S293>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_cont_n_T EnabledSubsystem_g;/* '<S292>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_cont_n_T EnabledSubsystem_p;/* '<S291>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_cont_p_T EnabledSubsystem_bk;/* '<S15>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_cont_p_T EnabledSubsystem_h;/* '<S14>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_b;/* '<S13>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_a;/* '<S12>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem;/* '<S11>/Enabled Subsystem' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_pid_control_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  ros_slros2_internal_block_Sub_T obj_p;/* '<S294>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_h;/* '<S293>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_h4;/* '<S292>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_hq;/* '<S291>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_n;/* '<S15>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_c;/* '<S14>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_m;/* '<S13>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_k;/* '<S12>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_f;/* '<S11>/SourceBlock' */
  real_T UnitDelay3_DSTATE;            /* '<Root>/Unit Delay3' */
  real_T UnitDelay4_DSTATE;            /* '<Root>/Unit Delay4' */
  real_T UnitDelay5_DSTATE;            /* '<Root>/Unit Delay5' */
  real_T UnitDelay6_DSTATE;            /* '<Root>/Unit Delay6' */
  real_T UnitDelay2_DSTATE;            /* '<Root>/Unit Delay2' */
  real_T Memory_PreviousInput[3];      /* '<S16>/Memory' */
  real_T Memory1_PreviousInput[3];     /* '<S16>/Memory1' */
  real_T NextOutput;                   /* '<S288>/White Noise' */
  real_T NextOutput_j[4];              /* '<S306>/White Noise' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<Root>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_g;               /* '<S16>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S16>/To Workspace1' */

  struct {
    void *LoggedData;
  } ToWorkspace10_PWORK;               /* '<S16>/To Workspace10' */

  struct {
    void *LoggedData;
  } ToWorkspace11_PWORK;               /* '<S16>/To Workspace11' */

  struct {
    void *LoggedData;
  } ToWorkspace12_PWORK;               /* '<S16>/To Workspace12' */

  struct {
    void *LoggedData;
  } ToWorkspace13_PWORK;               /* '<S16>/To Workspace13' */

  struct {
    void *LoggedData;
  } ToWorkspace14_PWORK;               /* '<S16>/To Workspace14' */

  struct {
    void *LoggedData;
  } ToWorkspace15_PWORK;               /* '<S16>/To Workspace15' */

  struct {
    void *LoggedData;
  } ToWorkspace16_PWORK;               /* '<S16>/To Workspace16' */

  struct {
    void *LoggedData;
  } ToWorkspace17_PWORK;               /* '<S16>/To Workspace17' */

  struct {
    void *LoggedData;
  } ToWorkspace18_PWORK;               /* '<S16>/To Workspace18' */

  struct {
    void *LoggedData;
  } ToWorkspace19_PWORK;               /* '<S16>/To Workspace19' */

  struct {
    void *LoggedData;
  } ToWorkspace2_PWORK;                /* '<S16>/To Workspace2' */

  struct {
    void *LoggedData;
  } ToWorkspace20_PWORK;               /* '<S16>/To Workspace20' */

  struct {
    void *LoggedData;
  } ToWorkspace21_PWORK;               /* '<S16>/To Workspace21' */

  struct {
    void *LoggedData;
  } ToWorkspace22_PWORK;               /* '<S16>/To Workspace22' */

  struct {
    void *LoggedData;
  } ToWorkspace23_PWORK;               /* '<S16>/To Workspace23' */

  struct {
    void *LoggedData;
  } ToWorkspace24_PWORK;               /* '<S16>/To Workspace24' */

  struct {
    void *LoggedData;
  } ToWorkspace25_PWORK;               /* '<S16>/To Workspace25' */

  struct {
    void *LoggedData;
  } ToWorkspace26_PWORK;               /* '<S16>/To Workspace26' */

  struct {
    void *LoggedData;
  } ToWorkspace27_PWORK;               /* '<S16>/To Workspace27' */

  struct {
    void *LoggedData;
  } ToWorkspace28_PWORK;               /* '<S16>/To Workspace28' */

  struct {
    void *LoggedData;
  } ToWorkspace29_PWORK;               /* '<S16>/To Workspace29' */

  struct {
    void *LoggedData;
  } ToWorkspace3_PWORK;                /* '<S16>/To Workspace3' */

  struct {
    void *LoggedData;
  } ToWorkspace30_PWORK;               /* '<S16>/To Workspace30' */

  struct {
    void *LoggedData;
  } ToWorkspace31_PWORK;               /* '<S16>/To Workspace31' */

  struct {
    void *LoggedData;
  } ToWorkspace32_PWORK;               /* '<S16>/To Workspace32' */

  struct {
    void *LoggedData;
  } ToWorkspace33_PWORK;               /* '<S16>/To Workspace33' */

  struct {
    void *LoggedData;
  } ToWorkspace4_PWORK;                /* '<S16>/To Workspace4' */

  struct {
    void *LoggedData;
  } ToWorkspace5_PWORK;                /* '<S16>/To Workspace5' */

  struct {
    void *LoggedData;
  } ToWorkspace6_PWORK;                /* '<S16>/To Workspace6' */

  struct {
    void *LoggedData;
  } ToWorkspace7_PWORK;                /* '<S16>/To Workspace7' */

  struct {
    void *LoggedData;
  } ToWorkspace8_PWORK;                /* '<S16>/To Workspace8' */

  struct {
    void *LoggedData;
  } ToWorkspace9_PWORK;                /* '<S16>/To Workspace9' */

  uint32_T PreLookUpIndexSearchprobofexcee;
                        /* '<S313>/PreLook-Up Index Search  (prob of exceed)' */
  uint32_T RandSeed;                   /* '<S288>/White Noise' */
  uint32_T PreLookUpIndexSearchaltitude_DW;
                              /* '<S313>/PreLook-Up Index Search  (altitude)' */
  uint32_T RandSeed_i[4];              /* '<S306>/White Noise' */
  robotics_slcore_internal_bloc_T obj_cl;
                             /* '<Root>/Coordinate Transformation Conversion' */
  int8_T ifHeightMaxlowaltitudeelseifHei;
  /* '<S301>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  int8_T ifHeightMaxlowaltitudeelseifH_k;
  /* '<S302>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  boolean_T Memory_PreviousInput_o;    /* '<S43>/Memory' */
  boolean_T Memory_PreviousInput_m;    /* '<S97>/Memory' */
  boolean_T Memory_PreviousInput_a;    /* '<S255>/Memory' */
  boolean_T objisempty;                /* '<S294>/SourceBlock' */
  boolean_T objisempty_l;              /* '<S293>/SourceBlock' */
  boolean_T objisempty_c;              /* '<S292>/SourceBlock' */
  boolean_T objisempty_a;              /* '<S291>/SourceBlock' */
  boolean_T objisempty_k;              /* '<S15>/SourceBlock' */
  boolean_T objisempty_g;              /* '<S14>/SourceBlock' */
  boolean_T objisempty_gq;             /* '<S13>/SourceBlock' */
  boolean_T objisempty_g5;             /* '<S12>/SourceBlock' */
  boolean_T objisempty_f;              /* '<S11>/SourceBlock' */
  boolean_T objisempty_d;    /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_ft;             /* '<S2>/ServiceCaller' */
  boolean_T Hwgws_MODE;                /* '<S297>/Hwgw(s)' */
  boolean_T Hvgws_MODE;                /* '<S297>/Hvgw(s)' */
  boolean_T Hugws_MODE;                /* '<S297>/Hugw(s)' */
  boolean_T Hrgw_MODE;                 /* '<S296>/Hrgw' */
  boolean_T Hqgw_MODE;                 /* '<S296>/Hqgw' */
  boolean_T Hpgw_MODE;                 /* '<S296>/Hpgw' */
};

/* Continuous states (default storage) */
struct X_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S16>/Integrator' */
  real_T Integrator_CSTATE_n;          /* '<S107>/Integrator' */
  real_T Filter_CSTATE;                /* '<S102>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S53>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S48>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S159>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S154>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S211>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S206>/Filter' */
  real_T Integrator_CSTATE_f;          /* '<S265>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S260>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S16>/Integrator1' */
  real_T TransferFcn_CSTATE[2];        /* '<S16>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S16>/Transfer Fcn1' */
  real_T wg_p1_CSTATE[2];              /* '<S312>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S312>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S311>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S311>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S310>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S309>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S308>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S307>/pgw_p' */
};

/* State derivatives (default storage) */
struct XDot_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S16>/Integrator' */
  real_T Integrator_CSTATE_n;          /* '<S107>/Integrator' */
  real_T Filter_CSTATE;                /* '<S102>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S53>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S48>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S159>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S154>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S211>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S206>/Filter' */
  real_T Integrator_CSTATE_f;          /* '<S265>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S260>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S16>/Integrator1' */
  real_T TransferFcn_CSTATE[2];        /* '<S16>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S16>/Transfer Fcn1' */
  real_T wg_p1_CSTATE[2];              /* '<S312>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S312>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S311>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S311>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S310>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S309>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S308>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S307>/pgw_p' */
};

/* State disabled  */
struct XDis_pid_control_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S16>/Integrator' */
  boolean_T Integrator_CSTATE_n;       /* '<S107>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S102>/Filter' */
  boolean_T Integrator_CSTATE_m;       /* '<S53>/Integrator' */
  boolean_T Filter_CSTATE_g;           /* '<S48>/Filter' */
  boolean_T Integrator_CSTATE_p;       /* '<S159>/Integrator' */
  boolean_T Filter_CSTATE_m;           /* '<S154>/Filter' */
  boolean_T Integrator_CSTATE_d;       /* '<S211>/Integrator' */
  boolean_T Filter_CSTATE_f;           /* '<S206>/Filter' */
  boolean_T Integrator_CSTATE_f;       /* '<S265>/Integrator' */
  boolean_T Filter_CSTATE_l;           /* '<S260>/Filter' */
  boolean_T Integrator1_CSTATE;        /* '<S16>/Integrator1' */
  boolean_T TransferFcn_CSTATE[2];     /* '<S16>/Transfer Fcn' */
  boolean_T TransferFcn1_CSTATE;       /* '<S16>/Transfer Fcn1' */
  boolean_T wg_p1_CSTATE[2];           /* '<S312>/wg_p1' */
  boolean_T wg_p2_CSTATE[2];           /* '<S312>/wg_p2' */
  boolean_T vg_p1_CSTATE[2];           /* '<S311>/vg_p1' */
  boolean_T vgw_p2_CSTATE[2];          /* '<S311>/vgw_p2' */
  boolean_T ug_p_CSTATE[2];            /* '<S310>/ug_p' */
  boolean_T rgw_p_CSTATE[2];           /* '<S309>/rgw_p' */
  boolean_T qgw_p_CSTATE[2];           /* '<S308>/qgw_p' */
  boolean_T pgw_p_CSTATE[2];           /* '<S307>/pgw_p' */
};

/* Invariant block signals (default storage) */
struct ConstB_pid_control_V1_T {
  real_T UnitConversion;               /* '<S295>/Unit Conversion' */
  real_T UnitConversion_k;             /* '<S305>/Unit Conversion' */
  real_T sigma_wg;                     /* '<S314>/sigma_wg ' */
  real_T UnitConversion_n;             /* '<S299>/Unit Conversion' */
  real_T UnitConversion_c;             /* '<S333>/Unit Conversion' */
  real_T PreLookUpIndexSearchprobofe;
                        /* '<S313>/PreLook-Up Index Search  (prob of exceed)' */
  real_T Sqrt[4];                      /* '<S306>/Sqrt' */
  real_T Sqrt1;                        /* '<S306>/Sqrt1' */
  real_T Divide[4];                    /* '<S306>/Divide' */
  real_T motorspeed;                   /* '<S16>/Gain2' */
  real_T Sum;                          /* '<S323>/Sum' */
  real_T Sum_a;                        /* '<S315>/Sum' */
  real_T sqrt_a;                       /* '<S312>/sqrt' */
  real_T w4;                           /* '<S307>/w4' */
  real_T u16;                          /* '<S307>/u^1//6' */
  uint32_T PreLookUpIndexSearchprobo_g;
                        /* '<S313>/PreLook-Up Index Search  (prob of exceed)' */
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
   * Referenced by: '<S16>/Integrator'
   */
  real_T Integrator_IC[12];

  /* Expression: h_vec
   * Referenced by: '<S313>/PreLook-Up Index Search  (altitude)'
   */
  real_T PreLookUpIndexSearchaltitude_Br[12];

  /* Expression: sigma_vec'
   * Referenced by: '<S313>/Medium//High Altitude Intensity'
   */
  real_T MediumHighAltitudeIntensity_Tab[84];

  /* Computed Parameter: MediumHighAltitudeIntensity_max
   * Referenced by: '<S313>/Medium//High Altitude Intensity'
   */
  uint32_T MediumHighAltitudeIntensity_max[2];
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
  real_T odeY[42];
  real_T odeF[3][42];
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

  /* private member function(s) for subsystem '<S11>/Enabled Subsystem'*/
  static void pid_contr_EnabledSubsystem_Init(B_EnabledSubsystem_pid_contro_T
    *localB);
  static void pid_control_V1_EnabledSubsystem(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Float64 *rtu_In1, B_EnabledSubsystem_pid_contro_T *localB);

  /* private member function(s) for subsystem '<S14>/Enabled Subsystem'*/
  static void pid_con_EnabledSubsystem_d_Init(B_EnabledSubsystem_pid_cont_p_T
    *localB);
  static void pid_control__EnabledSubsystem_h(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Int64 *rtu_In1, B_EnabledSubsystem_pid_cont_p_T *localB);

  /* private member function(s) for subsystem '<S291>/Enabled Subsystem'*/
  static void pid_con_EnabledSubsystem_i_Init(B_EnabledSubsystem_pid_cont_n_T
    *localB);
  static void pid_control__EnabledSubsystem_p(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Bool *rtu_In1, B_EnabledSubsystem_pid_cont_n_T *localB);

  /* private member function(s) for subsystem '<Root>'*/
  void pid_con_Subscriber_setupImpl_on(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_contro_Subscriber_setupImpl(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_co_Subscriber_setupImpl_onh(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_c_Subscriber_setupImpl_onhg(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_cont_Subscriber_setupImpl_o(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_con_ServiceCaller_setupImpl(const ros_slros2_internal_block_Ser_T
    *obj);
  void pid__Subscriber_setupImpl_onhgd(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_Subscriber_setupImpl_onhgd0(const ros_slros2_internal_block_Sub_T
    *obj);
  void pi_Subscriber_setupImpl_onhgd03(const ros_slros2_internal_block_Sub_T
    *obj);
  void p_Subscriber_setupImpl_onhgd03r(const ros_slros2_internal_block_Sub_T
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
 * Block '<S112>/Proportional Gain' : Eliminated nontunable gain of 1
 * Block '<S151>/Kb' : Eliminated nontunable gain of 1
 * Block '<S289>/Cast' : Eliminate redundant data type conversion
 * Block '<S289>/Cast To Double' : Eliminate redundant data type conversion
 * Block '<S289>/Cast To Double1' : Eliminate redundant data type conversion
 * Block '<S289>/Cast To Double2' : Eliminate redundant data type conversion
 * Block '<S289>/Cast To Double3' : Eliminate redundant data type conversion
 * Block '<S289>/Cast To Double4' : Eliminate redundant data type conversion
 * Block '<S319>/Reshape' : Reshape block reduction
 * Block '<S319>/Reshape1' : Reshape block reduction
 * Block '<S321>/Reshape' : Reshape block reduction
 * Block '<S327>/Reshape' : Reshape block reduction
 * Block '<S327>/Reshape1' : Reshape block reduction
 * Block '<S329>/Reshape' : Reshape block reduction
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
 * '<S5>'   : 'pid_control_V1/MATLAB Function2'
 * '<S6>'   : 'pid_control_V1/PID ALERON'
 * '<S7>'   : 'pid_control_V1/PID ALTURA'
 * '<S8>'   : 'pid_control_V1/PID PITCH//ELEVATOR'
 * '<S9>'   : 'pid_control_V1/PID TIIMON'
 * '<S10>'  : 'pid_control_V1/PID VELOCIDAD'
 * '<S11>'  : 'pid_control_V1/Subscribe'
 * '<S12>'  : 'pid_control_V1/Subscribe-YAW'
 * '<S13>'  : 'pid_control_V1/Subscribe-YAW1'
 * '<S14>'  : 'pid_control_V1/Subscribe1'
 * '<S15>'  : 'pid_control_V1/Subscribe2'
 * '<S16>'  : 'pid_control_V1/Subsystem'
 * '<S17>'  : 'pid_control_V1/PID ALERON/Anti-windup'
 * '<S18>'  : 'pid_control_V1/PID ALERON/D Gain'
 * '<S19>'  : 'pid_control_V1/PID ALERON/External Derivative'
 * '<S20>'  : 'pid_control_V1/PID ALERON/Filter'
 * '<S21>'  : 'pid_control_V1/PID ALERON/Filter ICs'
 * '<S22>'  : 'pid_control_V1/PID ALERON/I Gain'
 * '<S23>'  : 'pid_control_V1/PID ALERON/Ideal P Gain'
 * '<S24>'  : 'pid_control_V1/PID ALERON/Ideal P Gain Fdbk'
 * '<S25>'  : 'pid_control_V1/PID ALERON/Integrator'
 * '<S26>'  : 'pid_control_V1/PID ALERON/Integrator ICs'
 * '<S27>'  : 'pid_control_V1/PID ALERON/N Copy'
 * '<S28>'  : 'pid_control_V1/PID ALERON/N Gain'
 * '<S29>'  : 'pid_control_V1/PID ALERON/P Copy'
 * '<S30>'  : 'pid_control_V1/PID ALERON/Parallel P Gain'
 * '<S31>'  : 'pid_control_V1/PID ALERON/Reset Signal'
 * '<S32>'  : 'pid_control_V1/PID ALERON/Saturation'
 * '<S33>'  : 'pid_control_V1/PID ALERON/Saturation Fdbk'
 * '<S34>'  : 'pid_control_V1/PID ALERON/Sum'
 * '<S35>'  : 'pid_control_V1/PID ALERON/Sum Fdbk'
 * '<S36>'  : 'pid_control_V1/PID ALERON/Tracking Mode'
 * '<S37>'  : 'pid_control_V1/PID ALERON/Tracking Mode Sum'
 * '<S38>'  : 'pid_control_V1/PID ALERON/Tsamp - Integral'
 * '<S39>'  : 'pid_control_V1/PID ALERON/Tsamp - Ngain'
 * '<S40>'  : 'pid_control_V1/PID ALERON/postSat Signal'
 * '<S41>'  : 'pid_control_V1/PID ALERON/preInt Signal'
 * '<S42>'  : 'pid_control_V1/PID ALERON/preSat Signal'
 * '<S43>'  : 'pid_control_V1/PID ALERON/Anti-windup/Cont. Clamping Parallel'
 * '<S44>'  : 'pid_control_V1/PID ALERON/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S45>'  : 'pid_control_V1/PID ALERON/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S46>'  : 'pid_control_V1/PID ALERON/D Gain/Internal Parameters'
 * '<S47>'  : 'pid_control_V1/PID ALERON/External Derivative/Error'
 * '<S48>'  : 'pid_control_V1/PID ALERON/Filter/Cont. Filter'
 * '<S49>'  : 'pid_control_V1/PID ALERON/Filter ICs/Internal IC - Filter'
 * '<S50>'  : 'pid_control_V1/PID ALERON/I Gain/Internal Parameters'
 * '<S51>'  : 'pid_control_V1/PID ALERON/Ideal P Gain/Passthrough'
 * '<S52>'  : 'pid_control_V1/PID ALERON/Ideal P Gain Fdbk/Disabled'
 * '<S53>'  : 'pid_control_V1/PID ALERON/Integrator/Continuous'
 * '<S54>'  : 'pid_control_V1/PID ALERON/Integrator ICs/Internal IC'
 * '<S55>'  : 'pid_control_V1/PID ALERON/N Copy/Disabled'
 * '<S56>'  : 'pid_control_V1/PID ALERON/N Gain/Internal Parameters'
 * '<S57>'  : 'pid_control_V1/PID ALERON/P Copy/Disabled'
 * '<S58>'  : 'pid_control_V1/PID ALERON/Parallel P Gain/Internal Parameters'
 * '<S59>'  : 'pid_control_V1/PID ALERON/Reset Signal/Disabled'
 * '<S60>'  : 'pid_control_V1/PID ALERON/Saturation/Enabled'
 * '<S61>'  : 'pid_control_V1/PID ALERON/Saturation Fdbk/Disabled'
 * '<S62>'  : 'pid_control_V1/PID ALERON/Sum/Sum_PID'
 * '<S63>'  : 'pid_control_V1/PID ALERON/Sum Fdbk/Disabled'
 * '<S64>'  : 'pid_control_V1/PID ALERON/Tracking Mode/Disabled'
 * '<S65>'  : 'pid_control_V1/PID ALERON/Tracking Mode Sum/Passthrough'
 * '<S66>'  : 'pid_control_V1/PID ALERON/Tsamp - Integral/TsSignalSpecification'
 * '<S67>'  : 'pid_control_V1/PID ALERON/Tsamp - Ngain/Passthrough'
 * '<S68>'  : 'pid_control_V1/PID ALERON/postSat Signal/Forward_Path'
 * '<S69>'  : 'pid_control_V1/PID ALERON/preInt Signal/Internal PreInt'
 * '<S70>'  : 'pid_control_V1/PID ALERON/preSat Signal/Forward_Path'
 * '<S71>'  : 'pid_control_V1/PID ALTURA/Anti-windup'
 * '<S72>'  : 'pid_control_V1/PID ALTURA/D Gain'
 * '<S73>'  : 'pid_control_V1/PID ALTURA/External Derivative'
 * '<S74>'  : 'pid_control_V1/PID ALTURA/Filter'
 * '<S75>'  : 'pid_control_V1/PID ALTURA/Filter ICs'
 * '<S76>'  : 'pid_control_V1/PID ALTURA/I Gain'
 * '<S77>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain'
 * '<S78>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain Fdbk'
 * '<S79>'  : 'pid_control_V1/PID ALTURA/Integrator'
 * '<S80>'  : 'pid_control_V1/PID ALTURA/Integrator ICs'
 * '<S81>'  : 'pid_control_V1/PID ALTURA/N Copy'
 * '<S82>'  : 'pid_control_V1/PID ALTURA/N Gain'
 * '<S83>'  : 'pid_control_V1/PID ALTURA/P Copy'
 * '<S84>'  : 'pid_control_V1/PID ALTURA/Parallel P Gain'
 * '<S85>'  : 'pid_control_V1/PID ALTURA/Reset Signal'
 * '<S86>'  : 'pid_control_V1/PID ALTURA/Saturation'
 * '<S87>'  : 'pid_control_V1/PID ALTURA/Saturation Fdbk'
 * '<S88>'  : 'pid_control_V1/PID ALTURA/Sum'
 * '<S89>'  : 'pid_control_V1/PID ALTURA/Sum Fdbk'
 * '<S90>'  : 'pid_control_V1/PID ALTURA/Tracking Mode'
 * '<S91>'  : 'pid_control_V1/PID ALTURA/Tracking Mode Sum'
 * '<S92>'  : 'pid_control_V1/PID ALTURA/Tsamp - Integral'
 * '<S93>'  : 'pid_control_V1/PID ALTURA/Tsamp - Ngain'
 * '<S94>'  : 'pid_control_V1/PID ALTURA/postSat Signal'
 * '<S95>'  : 'pid_control_V1/PID ALTURA/preInt Signal'
 * '<S96>'  : 'pid_control_V1/PID ALTURA/preSat Signal'
 * '<S97>'  : 'pid_control_V1/PID ALTURA/Anti-windup/Cont. Clamping Parallel'
 * '<S98>'  : 'pid_control_V1/PID ALTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S99>'  : 'pid_control_V1/PID ALTURA/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S100>' : 'pid_control_V1/PID ALTURA/D Gain/Internal Parameters'
 * '<S101>' : 'pid_control_V1/PID ALTURA/External Derivative/Error'
 * '<S102>' : 'pid_control_V1/PID ALTURA/Filter/Cont. Filter'
 * '<S103>' : 'pid_control_V1/PID ALTURA/Filter ICs/Internal IC - Filter'
 * '<S104>' : 'pid_control_V1/PID ALTURA/I Gain/Internal Parameters'
 * '<S105>' : 'pid_control_V1/PID ALTURA/Ideal P Gain/Passthrough'
 * '<S106>' : 'pid_control_V1/PID ALTURA/Ideal P Gain Fdbk/Disabled'
 * '<S107>' : 'pid_control_V1/PID ALTURA/Integrator/Continuous'
 * '<S108>' : 'pid_control_V1/PID ALTURA/Integrator ICs/Internal IC'
 * '<S109>' : 'pid_control_V1/PID ALTURA/N Copy/Disabled'
 * '<S110>' : 'pid_control_V1/PID ALTURA/N Gain/Internal Parameters'
 * '<S111>' : 'pid_control_V1/PID ALTURA/P Copy/Disabled'
 * '<S112>' : 'pid_control_V1/PID ALTURA/Parallel P Gain/Internal Parameters'
 * '<S113>' : 'pid_control_V1/PID ALTURA/Reset Signal/Disabled'
 * '<S114>' : 'pid_control_V1/PID ALTURA/Saturation/Enabled'
 * '<S115>' : 'pid_control_V1/PID ALTURA/Saturation Fdbk/Disabled'
 * '<S116>' : 'pid_control_V1/PID ALTURA/Sum/Sum_PID'
 * '<S117>' : 'pid_control_V1/PID ALTURA/Sum Fdbk/Disabled'
 * '<S118>' : 'pid_control_V1/PID ALTURA/Tracking Mode/Disabled'
 * '<S119>' : 'pid_control_V1/PID ALTURA/Tracking Mode Sum/Passthrough'
 * '<S120>' : 'pid_control_V1/PID ALTURA/Tsamp - Integral/TsSignalSpecification'
 * '<S121>' : 'pid_control_V1/PID ALTURA/Tsamp - Ngain/Passthrough'
 * '<S122>' : 'pid_control_V1/PID ALTURA/postSat Signal/Forward_Path'
 * '<S123>' : 'pid_control_V1/PID ALTURA/preInt Signal/Internal PreInt'
 * '<S124>' : 'pid_control_V1/PID ALTURA/preSat Signal/Forward_Path'
 * '<S125>' : 'pid_control_V1/PID PITCH//ELEVATOR/Anti-windup'
 * '<S126>' : 'pid_control_V1/PID PITCH//ELEVATOR/D Gain'
 * '<S127>' : 'pid_control_V1/PID PITCH//ELEVATOR/External Derivative'
 * '<S128>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter'
 * '<S129>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter ICs'
 * '<S130>' : 'pid_control_V1/PID PITCH//ELEVATOR/I Gain'
 * '<S131>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain'
 * '<S132>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain Fdbk'
 * '<S133>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator'
 * '<S134>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator ICs'
 * '<S135>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Copy'
 * '<S136>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Gain'
 * '<S137>' : 'pid_control_V1/PID PITCH//ELEVATOR/P Copy'
 * '<S138>' : 'pid_control_V1/PID PITCH//ELEVATOR/Parallel P Gain'
 * '<S139>' : 'pid_control_V1/PID PITCH//ELEVATOR/Reset Signal'
 * '<S140>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation'
 * '<S141>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation Fdbk'
 * '<S142>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum'
 * '<S143>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum Fdbk'
 * '<S144>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode'
 * '<S145>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode Sum'
 * '<S146>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Integral'
 * '<S147>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Ngain'
 * '<S148>' : 'pid_control_V1/PID PITCH//ELEVATOR/postSat Signal'
 * '<S149>' : 'pid_control_V1/PID PITCH//ELEVATOR/preInt Signal'
 * '<S150>' : 'pid_control_V1/PID PITCH//ELEVATOR/preSat Signal'
 * '<S151>' : 'pid_control_V1/PID PITCH//ELEVATOR/Anti-windup/Back Calculation'
 * '<S152>' : 'pid_control_V1/PID PITCH//ELEVATOR/D Gain/Internal Parameters'
 * '<S153>' : 'pid_control_V1/PID PITCH//ELEVATOR/External Derivative/Error'
 * '<S154>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter/Cont. Filter'
 * '<S155>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter ICs/Internal IC - Filter'
 * '<S156>' : 'pid_control_V1/PID PITCH//ELEVATOR/I Gain/Internal Parameters'
 * '<S157>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain/Passthrough'
 * '<S158>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain Fdbk/Disabled'
 * '<S159>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator/Continuous'
 * '<S160>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator ICs/Internal IC'
 * '<S161>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Copy/Disabled'
 * '<S162>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Gain/Internal Parameters'
 * '<S163>' : 'pid_control_V1/PID PITCH//ELEVATOR/P Copy/Disabled'
 * '<S164>' : 'pid_control_V1/PID PITCH//ELEVATOR/Parallel P Gain/Internal Parameters'
 * '<S165>' : 'pid_control_V1/PID PITCH//ELEVATOR/Reset Signal/Disabled'
 * '<S166>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation/Enabled'
 * '<S167>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation Fdbk/Disabled'
 * '<S168>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum/Sum_PID'
 * '<S169>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum Fdbk/Disabled'
 * '<S170>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode/Disabled'
 * '<S171>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode Sum/Passthrough'
 * '<S172>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Integral/TsSignalSpecification'
 * '<S173>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Ngain/Passthrough'
 * '<S174>' : 'pid_control_V1/PID PITCH//ELEVATOR/postSat Signal/Forward_Path'
 * '<S175>' : 'pid_control_V1/PID PITCH//ELEVATOR/preInt Signal/Internal PreInt'
 * '<S176>' : 'pid_control_V1/PID PITCH//ELEVATOR/preSat Signal/Forward_Path'
 * '<S177>' : 'pid_control_V1/PID TIIMON/Anti-windup'
 * '<S178>' : 'pid_control_V1/PID TIIMON/D Gain'
 * '<S179>' : 'pid_control_V1/PID TIIMON/External Derivative'
 * '<S180>' : 'pid_control_V1/PID TIIMON/Filter'
 * '<S181>' : 'pid_control_V1/PID TIIMON/Filter ICs'
 * '<S182>' : 'pid_control_V1/PID TIIMON/I Gain'
 * '<S183>' : 'pid_control_V1/PID TIIMON/Ideal P Gain'
 * '<S184>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk'
 * '<S185>' : 'pid_control_V1/PID TIIMON/Integrator'
 * '<S186>' : 'pid_control_V1/PID TIIMON/Integrator ICs'
 * '<S187>' : 'pid_control_V1/PID TIIMON/N Copy'
 * '<S188>' : 'pid_control_V1/PID TIIMON/N Gain'
 * '<S189>' : 'pid_control_V1/PID TIIMON/P Copy'
 * '<S190>' : 'pid_control_V1/PID TIIMON/Parallel P Gain'
 * '<S191>' : 'pid_control_V1/PID TIIMON/Reset Signal'
 * '<S192>' : 'pid_control_V1/PID TIIMON/Saturation'
 * '<S193>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk'
 * '<S194>' : 'pid_control_V1/PID TIIMON/Sum'
 * '<S195>' : 'pid_control_V1/PID TIIMON/Sum Fdbk'
 * '<S196>' : 'pid_control_V1/PID TIIMON/Tracking Mode'
 * '<S197>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum'
 * '<S198>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral'
 * '<S199>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain'
 * '<S200>' : 'pid_control_V1/PID TIIMON/postSat Signal'
 * '<S201>' : 'pid_control_V1/PID TIIMON/preInt Signal'
 * '<S202>' : 'pid_control_V1/PID TIIMON/preSat Signal'
 * '<S203>' : 'pid_control_V1/PID TIIMON/Anti-windup/Passthrough'
 * '<S204>' : 'pid_control_V1/PID TIIMON/D Gain/Internal Parameters'
 * '<S205>' : 'pid_control_V1/PID TIIMON/External Derivative/Error'
 * '<S206>' : 'pid_control_V1/PID TIIMON/Filter/Cont. Filter'
 * '<S207>' : 'pid_control_V1/PID TIIMON/Filter ICs/Internal IC - Filter'
 * '<S208>' : 'pid_control_V1/PID TIIMON/I Gain/Internal Parameters'
 * '<S209>' : 'pid_control_V1/PID TIIMON/Ideal P Gain/Passthrough'
 * '<S210>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk/Disabled'
 * '<S211>' : 'pid_control_V1/PID TIIMON/Integrator/Continuous'
 * '<S212>' : 'pid_control_V1/PID TIIMON/Integrator ICs/Internal IC'
 * '<S213>' : 'pid_control_V1/PID TIIMON/N Copy/Disabled'
 * '<S214>' : 'pid_control_V1/PID TIIMON/N Gain/Internal Parameters'
 * '<S215>' : 'pid_control_V1/PID TIIMON/P Copy/Disabled'
 * '<S216>' : 'pid_control_V1/PID TIIMON/Parallel P Gain/Internal Parameters'
 * '<S217>' : 'pid_control_V1/PID TIIMON/Reset Signal/Disabled'
 * '<S218>' : 'pid_control_V1/PID TIIMON/Saturation/Enabled'
 * '<S219>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk/Disabled'
 * '<S220>' : 'pid_control_V1/PID TIIMON/Sum/Sum_PID'
 * '<S221>' : 'pid_control_V1/PID TIIMON/Sum Fdbk/Disabled'
 * '<S222>' : 'pid_control_V1/PID TIIMON/Tracking Mode/Disabled'
 * '<S223>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum/Passthrough'
 * '<S224>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral/TsSignalSpecification'
 * '<S225>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain/Passthrough'
 * '<S226>' : 'pid_control_V1/PID TIIMON/postSat Signal/Forward_Path'
 * '<S227>' : 'pid_control_V1/PID TIIMON/preInt Signal/Internal PreInt'
 * '<S228>' : 'pid_control_V1/PID TIIMON/preSat Signal/Forward_Path'
 * '<S229>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup'
 * '<S230>' : 'pid_control_V1/PID VELOCIDAD/D Gain'
 * '<S231>' : 'pid_control_V1/PID VELOCIDAD/External Derivative'
 * '<S232>' : 'pid_control_V1/PID VELOCIDAD/Filter'
 * '<S233>' : 'pid_control_V1/PID VELOCIDAD/Filter ICs'
 * '<S234>' : 'pid_control_V1/PID VELOCIDAD/I Gain'
 * '<S235>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain'
 * '<S236>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain Fdbk'
 * '<S237>' : 'pid_control_V1/PID VELOCIDAD/Integrator'
 * '<S238>' : 'pid_control_V1/PID VELOCIDAD/Integrator ICs'
 * '<S239>' : 'pid_control_V1/PID VELOCIDAD/N Copy'
 * '<S240>' : 'pid_control_V1/PID VELOCIDAD/N Gain'
 * '<S241>' : 'pid_control_V1/PID VELOCIDAD/P Copy'
 * '<S242>' : 'pid_control_V1/PID VELOCIDAD/Parallel P Gain'
 * '<S243>' : 'pid_control_V1/PID VELOCIDAD/Reset Signal'
 * '<S244>' : 'pid_control_V1/PID VELOCIDAD/Saturation'
 * '<S245>' : 'pid_control_V1/PID VELOCIDAD/Saturation Fdbk'
 * '<S246>' : 'pid_control_V1/PID VELOCIDAD/Sum'
 * '<S247>' : 'pid_control_V1/PID VELOCIDAD/Sum Fdbk'
 * '<S248>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode'
 * '<S249>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode Sum'
 * '<S250>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Integral'
 * '<S251>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Ngain'
 * '<S252>' : 'pid_control_V1/PID VELOCIDAD/postSat Signal'
 * '<S253>' : 'pid_control_V1/PID VELOCIDAD/preInt Signal'
 * '<S254>' : 'pid_control_V1/PID VELOCIDAD/preSat Signal'
 * '<S255>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel'
 * '<S256>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S257>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S258>' : 'pid_control_V1/PID VELOCIDAD/D Gain/Internal Parameters'
 * '<S259>' : 'pid_control_V1/PID VELOCIDAD/External Derivative/Error'
 * '<S260>' : 'pid_control_V1/PID VELOCIDAD/Filter/Cont. Filter'
 * '<S261>' : 'pid_control_V1/PID VELOCIDAD/Filter ICs/Internal IC - Filter'
 * '<S262>' : 'pid_control_V1/PID VELOCIDAD/I Gain/Internal Parameters'
 * '<S263>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain/Passthrough'
 * '<S264>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain Fdbk/Disabled'
 * '<S265>' : 'pid_control_V1/PID VELOCIDAD/Integrator/Continuous'
 * '<S266>' : 'pid_control_V1/PID VELOCIDAD/Integrator ICs/Internal IC'
 * '<S267>' : 'pid_control_V1/PID VELOCIDAD/N Copy/Disabled'
 * '<S268>' : 'pid_control_V1/PID VELOCIDAD/N Gain/Internal Parameters'
 * '<S269>' : 'pid_control_V1/PID VELOCIDAD/P Copy/Disabled'
 * '<S270>' : 'pid_control_V1/PID VELOCIDAD/Parallel P Gain/Internal Parameters'
 * '<S271>' : 'pid_control_V1/PID VELOCIDAD/Reset Signal/Disabled'
 * '<S272>' : 'pid_control_V1/PID VELOCIDAD/Saturation/Enabled'
 * '<S273>' : 'pid_control_V1/PID VELOCIDAD/Saturation Fdbk/Disabled'
 * '<S274>' : 'pid_control_V1/PID VELOCIDAD/Sum/Sum_PID'
 * '<S275>' : 'pid_control_V1/PID VELOCIDAD/Sum Fdbk/Disabled'
 * '<S276>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode/Disabled'
 * '<S277>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode Sum/Passthrough'
 * '<S278>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Integral/TsSignalSpecification'
 * '<S279>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Ngain/Passthrough'
 * '<S280>' : 'pid_control_V1/PID VELOCIDAD/postSat Signal/Forward_Path'
 * '<S281>' : 'pid_control_V1/PID VELOCIDAD/preInt Signal/Internal PreInt'
 * '<S282>' : 'pid_control_V1/PID VELOCIDAD/preSat Signal/Forward_Path'
 * '<S283>' : 'pid_control_V1/Subscribe/Enabled Subsystem'
 * '<S284>' : 'pid_control_V1/Subscribe-YAW/Enabled Subsystem'
 * '<S285>' : 'pid_control_V1/Subscribe-YAW1/Enabled Subsystem'
 * '<S286>' : 'pid_control_V1/Subscribe1/Enabled Subsystem'
 * '<S287>' : 'pid_control_V1/Subscribe2/Enabled Subsystem'
 * '<S288>' : 'pid_control_V1/Subsystem/Band-Limited White Noise'
 * '<S289>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))'
 * '<S290>' : 'pid_control_V1/Subsystem/MATLAB Function'
 * '<S291>' : 'pid_control_V1/Subsystem/Subscribe'
 * '<S292>' : 'pid_control_V1/Subsystem/Subscribe1'
 * '<S293>' : 'pid_control_V1/Subsystem/Subscribe2'
 * '<S294>' : 'pid_control_V1/Subsystem/Subscribe3'
 * '<S295>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Angle Conversion'
 * '<S296>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates'
 * '<S297>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities'
 * '<S298>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Length Conversion'
 * '<S299>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Length Conversion1'
 * '<S300>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities'
 * '<S301>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates'
 * '<S302>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities'
 * '<S303>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths'
 * '<S304>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Velocity Conversion'
 * '<S305>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Velocity Conversion2'
 * '<S306>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/White Noise'
 * '<S307>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hpgw'
 * '<S308>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hqgw'
 * '<S309>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hrgw'
 * '<S310>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hugw(s)'
 * '<S311>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hvgw(s)'
 * '<S312>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hwgw(s)'
 * '<S313>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities/High Altitude Intensity'
 * '<S314>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities/Low Altitude Intensity'
 * '<S315>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates'
 * '<S316>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates'
 * '<S317>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Medium//High  altitude rates'
 * '<S318>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Merge Subsystems'
 * '<S319>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates/wind to body transformation'
 * '<S320>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates/wind to body transformation/convert to earth coords'
 * '<S321>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates/wind to body transformation'
 * '<S322>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates/wind to body transformation/convert to earth coords'
 * '<S323>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities'
 * '<S324>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities'
 * '<S325>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Medium//High  altitude velocities'
 * '<S326>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Merge Subsystems'
 * '<S327>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities/wind to body transformation'
 * '<S328>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities/wind to body transformation/convert to earth coords'
 * '<S329>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities/wind to body transformation'
 * '<S330>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities/wind to body transformation/convert to earth coords'
 * '<S331>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Low altitude scale length'
 * '<S332>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Medium//High altitude scale length'
 * '<S333>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Medium//High altitude scale length/Length Conversion'
 * '<S334>' : 'pid_control_V1/Subsystem/Subscribe/Enabled Subsystem'
 * '<S335>' : 'pid_control_V1/Subsystem/Subscribe1/Enabled Subsystem'
 * '<S336>' : 'pid_control_V1/Subsystem/Subscribe2/Enabled Subsystem'
 * '<S337>' : 'pid_control_V1/Subsystem/Subscribe3/Enabled Subsystem'
 */
#endif                                 /* pid_control_V1_h_ */
