/*
 * pid_control_V1.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "pid_control_V1".
 *
 * Model version              : 12.59
 * Simulink Coder version : 25.2 (R2025b) 28-Jul-2025
 * C++ source code generated on : Mon Mar  9 19:59:35 2026
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

/* Block signals for system '<S10>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_contro_T {
  SL_Bus_std_msgs_Float64 In1;         /* '<S275>/In1' */
};

/* Block signals for system '<S279>/Enabled Subsystem' */
struct B_EnabledSubsystem_pid_cont_n_T {
  SL_Bus_std_msgs_Bool In1;            /* '<S320>/In1' */
};

/* Block signals (default storage) */
struct B_pid_control_V1_T {
  SL_Bus_gazebo_msgs_SetEntityStateRequest BusAssignment;/* '<Root>/Bus Assignment' */
  real_T x[12];                        /* '<S12>/Integrator' */
  uint8_T stringOut[128];              /* '<Root>/MATLAB Function1' */
  uint8_T stringOut_l[128];            /* '<Root>/MATLAB Function' */
  real_T RotationAnglestoDirectionCo[9];
                        /* '<S12>/Rotation Angles to Direction Cosine Matrix' */
  real_T FA_b_tmp_p[9];
  real_T TmpSignalConversionAtSFunct[5];/* '<S12>/MATLAB Function' */
  char_T b_zeroDelimTopic[25];
  real_T Memory[3];                    /* '<S12>/Memory' */
  real_T Memory1[3];                   /* '<S12>/Memory1' */
  real_T wbe_b[3];
  real_T FA_b[3];
  real_T Product_be[3];                /* '<S305>/Product' */
  real_T Dtot[3];
  char_T b_zeroDelimTopic_m[22];
  char_T b_zeroDelimTopic_c[22];
  char_T b_zeroDelimTopic_k[17];
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  real_T frac[2];
  real_T dv[2];
  real_T Switch3;                      /* '<Root>/Switch3' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T FilterCoefficient;            /* '<S102>/Filter Coefficient' */
  real_T Sum;                          /* '<S108>/Sum' */
  real_T Saturation1;                  /* '<Root>/Saturation1' */
  real_T FilterCoefficient_c;          /* '<S50>/Filter Coefficient' */
  real_T Saturation;                   /* '<S54>/Saturation' */
  real_T Switch2;                      /* '<Root>/Switch2' */
  real_T FilterCoefficient_p;          /* '<S206>/Filter Coefficient' */
  real_T Saturation_m;                 /* '<S210>/Saturation' */
  real_T Saturation_i;                 /* '<Root>/Saturation' */
  real_T FilterCoefficient_m;          /* '<S154>/Filter Coefficient' */
  real_T Saturation_f;                 /* '<S158>/Saturation' */
  real_T SumI4;                        /* '<S39>/SumI4' */
  real_T IntegralGain;                 /* '<S96>/Integral Gain' */
  real_T SumI4_i;                      /* '<S143>/SumI4' */
  real_T IntegralGain_a;               /* '<S200>/Integral Gain' */
  real_T FilterCoefficient_cv;         /* '<S260>/Filter Coefficient' */
  real_T Switch;                       /* '<S247>/Switch' */
  real_T Saturation_o;                 /* '<S264>/Saturation' */
  real_T Power;                        /* '<S12>/Product2' */
  real_T Gain3;                        /* '<S12>/Gain3' */
  real_T EnergykWh;                    /* '<S12>/Gain1' */
  real_T powerdemand;                  /* '<S12>/Divide' */
  real_T loadtorque;                   /* '<S12>/Divide1' */
  real_T Product[4];                   /* '<S292>/Product' */
  real_T Switch_p[3];                  /* '<S12>/Switch' */
  real_T Switch1[3];                   /* '<S12>/Switch1' */
  real_T XDOT[40];                     /* '<S12>/MATLAB Function' */
  real_T w[2];                         /* '<S298>/w' */
  real_T w_a[2];                       /* '<S298>/w ' */
  real_T LwgV1[2];                     /* '<S298>/Lwg//V 1' */
  real_T w_g[2];                       /* '<S297>/w' */
  real_T w_e[2];                       /* '<S297>/w ' */
  real_T w1[2];                        /* '<S297>/w 1' */
  real_T w_n[2];                       /* '<S296>/w' */
  real_T w1_c[2];                      /* '<S296>/w1' */
  real_T w_d[2];                       /* '<S295>/w' */
  real_T w_e0[2];                      /* '<S294>/w' */
  real_T UnaryMinus[2];                /* '<S294>/Unary Minus' */
  real_T w_o[2];                       /* '<S293>/w' */
  real_T sigma_w[2];                   /* '<S293>/sigma_w' */
  real_T u2;
  real_T Va;
  real_T Q;
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
  real_T Sum2_l;                       /* '<Root>/Sum2' */
  real_T Sum_ks;                       /* '<S56>/Sum' */
  real_T Sum4;                         /* '<Root>/Sum4' */
  real_T Sum5;                         /* '<Root>/Sum5' */
  real_T Sum_hl;                       /* '<S160>/Sum' */
  real_T Sum1_g;                       /* '<Root>/Sum1' */
  real_T FA_b_tmp;
  real_T FE1_b_idx_1;
  real_T Mcg_b_idx_2;
  real_T Mcg_b_idx_0;
  real_T FE2_b_idx_0;
  real_T FE2_b_idx_2;
  real_T Fg_b_idx_2;
  real_T Fg_b_idx_1;
  real_T FA_b_c;
  real_T c_the_tmp;
  real_T c_the_tmp_b;
  uint32_T bpIndex[2];
  uint32_T lengthOut;                  /* '<Root>/MATLAB Function1' */
  uint32_T lengthOut_e;                /* '<Root>/MATLAB Function' */
  boolean_T AND3;                      /* '<S247>/AND3' */
  boolean_T Memory_h;                  /* '<S247>/Memory' */
  boolean_T SourceBlock_o1;            /* '<S280>/SourceBlock' */
  boolean_T SourceBlock_o1_h;          /* '<S279>/SourceBlock' */
  boolean_T SourceBlock_o1_o;          /* '<S11>/SourceBlock' */
  boolean_T SourceBlock_o1_g;          /* '<S10>/SourceBlock' */
  B_EnabledSubsystem_pid_cont_n_T EnabledSubsystem_g;/* '<S280>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_cont_n_T EnabledSubsystem_p;/* '<S279>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem_b;/* '<S11>/Enabled Subsystem' */
  B_EnabledSubsystem_pid_contro_T EnabledSubsystem;/* '<S10>/Enabled Subsystem' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_pid_control_V1_T {
  ros_slros2_internal_block_Ser_T obj; /* '<S2>/ServiceCaller' */
  ros_slros2_internal_block_Sub_T obj_h;/* '<S280>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_hq;/* '<S279>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_m;/* '<S11>/SourceBlock' */
  ros_slros2_internal_block_Sub_T obj_k;/* '<S10>/SourceBlock' */
  real_T UnitDelay3_DSTATE;            /* '<Root>/Unit Delay3' */
  real_T UnitDelay2_DSTATE;            /* '<Root>/Unit Delay2' */
  real_T Memory_PreviousInput[3];      /* '<S12>/Memory' */
  real_T Memory1_PreviousInput[3];     /* '<S12>/Memory1' */
  real_T NextOutput[4];                /* '<S292>/White Noise' */
  struct {
    void *LoggedData;
  } ToWorkspace_PWORK;                 /* '<Root>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace_PWORK_g;               /* '<S12>/To Workspace' */

  struct {
    void *LoggedData;
  } ToWorkspace1_PWORK;                /* '<S12>/To Workspace1' */

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
  } ToWorkspace30_PWORK;               /* '<S12>/To Workspace30' */

  struct {
    void *LoggedData;
  } ToWorkspace31_PWORK;               /* '<S12>/To Workspace31' */

  struct {
    void *LoggedData;
  } ToWorkspace32_PWORK;               /* '<S12>/To Workspace32' */

  struct {
    void *LoggedData;
  } ToWorkspace33_PWORK;               /* '<S12>/To Workspace33' */

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

  uint32_T PreLookUpIndexSearchprobofexcee;
                        /* '<S299>/PreLook-Up Index Search  (prob of exceed)' */
  uint32_T PreLookUpIndexSearchaltitude_DW;
                              /* '<S299>/PreLook-Up Index Search  (altitude)' */
  uint32_T RandSeed[4];                /* '<S292>/White Noise' */
  robotics_slcore_internal_bloc_T obj_c;
                             /* '<Root>/Coordinate Transformation Conversion' */
  int8_T ifHeightMaxlowaltitudeelseifHei;
  /* '<S287>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  int8_T ifHeightMaxlowaltitudeelseifH_k;
  /* '<S288>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  boolean_T Memory_PreviousInput_a;    /* '<S247>/Memory' */
  boolean_T objisempty;                /* '<S280>/SourceBlock' */
  boolean_T objisempty_a;              /* '<S279>/SourceBlock' */
  boolean_T objisempty_g;              /* '<S11>/SourceBlock' */
  boolean_T objisempty_g5;             /* '<S10>/SourceBlock' */
  boolean_T objisempty_d;    /* '<Root>/Coordinate Transformation Conversion' */
  boolean_T objisempty_f;              /* '<S2>/ServiceCaller' */
  boolean_T Hwgws_MODE;                /* '<S283>/Hwgw(s)' */
  boolean_T Hvgws_MODE;                /* '<S283>/Hvgw(s)' */
  boolean_T Hugws_MODE;                /* '<S283>/Hugw(s)' */
  boolean_T Hrgw_MODE;                 /* '<S282>/Hrgw' */
  boolean_T Hqgw_MODE;                 /* '<S282>/Hqgw' */
  boolean_T Hpgw_MODE;                 /* '<S282>/Hpgw' */
};

/* Continuous states (default storage) */
struct X_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_n;          /* '<S99>/Integrator' */
  real_T Filter_CSTATE;                /* '<S94>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S47>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S42>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S203>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S198>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S151>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S146>/Filter' */
  real_T Integrator_CSTATE_f;          /* '<S257>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S252>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S12>/Integrator1' */
  real_T wg_p1_CSTATE[2];              /* '<S298>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S298>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S297>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S297>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S296>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S295>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S294>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S293>/pgw_p' */
};

/* State derivatives (default storage) */
struct XDot_pid_control_V1_T {
  real_T Integrator_CSTATE[12];        /* '<S12>/Integrator' */
  real_T Integrator_CSTATE_n;          /* '<S99>/Integrator' */
  real_T Filter_CSTATE;                /* '<S94>/Filter' */
  real_T Integrator_CSTATE_m;          /* '<S47>/Integrator' */
  real_T Filter_CSTATE_g;              /* '<S42>/Filter' */
  real_T Integrator_CSTATE_d;          /* '<S203>/Integrator' */
  real_T Filter_CSTATE_f;              /* '<S198>/Filter' */
  real_T Integrator_CSTATE_p;          /* '<S151>/Integrator' */
  real_T Filter_CSTATE_m;              /* '<S146>/Filter' */
  real_T Integrator_CSTATE_f;          /* '<S257>/Integrator' */
  real_T Filter_CSTATE_l;              /* '<S252>/Filter' */
  real_T Integrator1_CSTATE;           /* '<S12>/Integrator1' */
  real_T wg_p1_CSTATE[2];              /* '<S298>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S298>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S297>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S297>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S296>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S295>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S294>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S293>/pgw_p' */
};

/* State disabled  */
struct XDis_pid_control_V1_T {
  boolean_T Integrator_CSTATE[12];     /* '<S12>/Integrator' */
  boolean_T Integrator_CSTATE_n;       /* '<S99>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S94>/Filter' */
  boolean_T Integrator_CSTATE_m;       /* '<S47>/Integrator' */
  boolean_T Filter_CSTATE_g;           /* '<S42>/Filter' */
  boolean_T Integrator_CSTATE_d;       /* '<S203>/Integrator' */
  boolean_T Filter_CSTATE_f;           /* '<S198>/Filter' */
  boolean_T Integrator_CSTATE_p;       /* '<S151>/Integrator' */
  boolean_T Filter_CSTATE_m;           /* '<S146>/Filter' */
  boolean_T Integrator_CSTATE_f;       /* '<S257>/Integrator' */
  boolean_T Filter_CSTATE_l;           /* '<S252>/Filter' */
  boolean_T Integrator1_CSTATE;        /* '<S12>/Integrator1' */
  boolean_T wg_p1_CSTATE[2];           /* '<S298>/wg_p1' */
  boolean_T wg_p2_CSTATE[2];           /* '<S298>/wg_p2' */
  boolean_T vg_p1_CSTATE[2];           /* '<S297>/vg_p1' */
  boolean_T vgw_p2_CSTATE[2];          /* '<S297>/vgw_p2' */
  boolean_T ug_p_CSTATE[2];            /* '<S296>/ug_p' */
  boolean_T rgw_p_CSTATE[2];           /* '<S295>/rgw_p' */
  boolean_T qgw_p_CSTATE[2];           /* '<S294>/qgw_p' */
  boolean_T pgw_p_CSTATE[2];           /* '<S293>/pgw_p' */
};

/* Invariant block signals (default storage) */
struct ConstB_pid_control_V1_T {
  real_T UnitConversion;               /* '<S281>/Unit Conversion' */
  real_T UnitConversion_k;             /* '<S291>/Unit Conversion' */
  real_T sigma_wg;                     /* '<S300>/sigma_wg ' */
  real_T UnitConversion_n;             /* '<S285>/Unit Conversion' */
  real_T UnitConversion_c;             /* '<S319>/Unit Conversion' */
  real_T PreLookUpIndexSearchprobofe;
                        /* '<S299>/PreLook-Up Index Search  (prob of exceed)' */
  real_T Sqrt[4];                      /* '<S292>/Sqrt' */
  real_T Sqrt1;                        /* '<S292>/Sqrt1' */
  real_T Divide[4];                    /* '<S292>/Divide' */
  real_T motorspeed;                   /* '<S12>/Gain2' */
  real_T Sum;                          /* '<S309>/Sum' */
  real_T Sum_a;                        /* '<S301>/Sum' */
  real_T sqrt_a;                       /* '<S298>/sqrt' */
  real_T w4;                           /* '<S293>/w4' */
  real_T u16;                          /* '<S293>/u^1//6' */
  uint32_T PreLookUpIndexSearchprobo_g;
                        /* '<S299>/PreLook-Up Index Search  (prob of exceed)' */
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

  /* Expression: h_vec
   * Referenced by: '<S299>/PreLook-Up Index Search  (altitude)'
   */
  real_T PreLookUpIndexSearchaltitude_Br[12];

  /* Expression: sigma_vec'
   * Referenced by: '<S299>/Medium//High Altitude Intensity'
   */
  real_T MediumHighAltitudeIntensity_Tab[84];

  /* Computed Parameter: MediumHighAltitudeIntensity_max
   * Referenced by: '<S299>/Medium//High Altitude Intensity'
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
  real_T odeY[39];
  real_T odeF[3][39];
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

  /* private member function(s) for subsystem '<S10>/Enabled Subsystem'*/
  static void pid_contr_EnabledSubsystem_Init(B_EnabledSubsystem_pid_contro_T
    *localB);
  static void pid_control_V1_EnabledSubsystem(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Float64 *rtu_In1, B_EnabledSubsystem_pid_contro_T *localB);

  /* private member function(s) for subsystem '<S279>/Enabled Subsystem'*/
  static void pid_con_EnabledSubsystem_i_Init(B_EnabledSubsystem_pid_cont_n_T
    *localB);
  static void pid_control__EnabledSubsystem_p(boolean_T rtu_Enable, const
    SL_Bus_std_msgs_Bool *rtu_In1, B_EnabledSubsystem_pid_cont_n_T *localB);

  /* private member function(s) for subsystem '<Root>'*/
  void pid_cont_Subscriber_setupImpl_o(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_contro_Subscriber_setupImpl(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_con_ServiceCaller_setupImpl(const ros_slros2_internal_block_Ser_T
    *obj);
  void pid_con_Subscriber_setupImpl_on(const ros_slros2_internal_block_Sub_T
    *obj);
  void pid_co_Subscriber_setupImpl_onh(const ros_slros2_internal_block_Sub_T
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
 * Block '<S39>/Kb' : Eliminated nontunable gain of 1
 * Block '<S104>/Proportional Gain' : Eliminated nontunable gain of 1
 * Block '<S143>/Kb' : Eliminated nontunable gain of 1
 * Block '<S277>/Cast' : Eliminate redundant data type conversion
 * Block '<S277>/Cast To Double' : Eliminate redundant data type conversion
 * Block '<S277>/Cast To Double1' : Eliminate redundant data type conversion
 * Block '<S277>/Cast To Double2' : Eliminate redundant data type conversion
 * Block '<S277>/Cast To Double3' : Eliminate redundant data type conversion
 * Block '<S277>/Cast To Double4' : Eliminate redundant data type conversion
 * Block '<S305>/Reshape' : Reshape block reduction
 * Block '<S305>/Reshape1' : Reshape block reduction
 * Block '<S307>/Reshape' : Reshape block reduction
 * Block '<S313>/Reshape' : Reshape block reduction
 * Block '<S313>/Reshape1' : Reshape block reduction
 * Block '<S315>/Reshape' : Reshape block reduction
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
 * '<S6>'   : 'pid_control_V1/PID ALTURA'
 * '<S7>'   : 'pid_control_V1/PID PITCH//ELEVATOR'
 * '<S8>'   : 'pid_control_V1/PID TIIMON'
 * '<S9>'   : 'pid_control_V1/PID VELOCIDAD'
 * '<S10>'  : 'pid_control_V1/Subscribe-YAW'
 * '<S11>'  : 'pid_control_V1/Subscribe-YAW1'
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
 * '<S39>'  : 'pid_control_V1/PID ALERON/Anti-windup/Back Calculation'
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
 * '<S65>'  : 'pid_control_V1/PID ALTURA/Anti-windup'
 * '<S66>'  : 'pid_control_V1/PID ALTURA/D Gain'
 * '<S67>'  : 'pid_control_V1/PID ALTURA/External Derivative'
 * '<S68>'  : 'pid_control_V1/PID ALTURA/Filter'
 * '<S69>'  : 'pid_control_V1/PID ALTURA/Filter ICs'
 * '<S70>'  : 'pid_control_V1/PID ALTURA/I Gain'
 * '<S71>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain'
 * '<S72>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain Fdbk'
 * '<S73>'  : 'pid_control_V1/PID ALTURA/Integrator'
 * '<S74>'  : 'pid_control_V1/PID ALTURA/Integrator ICs'
 * '<S75>'  : 'pid_control_V1/PID ALTURA/N Copy'
 * '<S76>'  : 'pid_control_V1/PID ALTURA/N Gain'
 * '<S77>'  : 'pid_control_V1/PID ALTURA/P Copy'
 * '<S78>'  : 'pid_control_V1/PID ALTURA/Parallel P Gain'
 * '<S79>'  : 'pid_control_V1/PID ALTURA/Reset Signal'
 * '<S80>'  : 'pid_control_V1/PID ALTURA/Saturation'
 * '<S81>'  : 'pid_control_V1/PID ALTURA/Saturation Fdbk'
 * '<S82>'  : 'pid_control_V1/PID ALTURA/Sum'
 * '<S83>'  : 'pid_control_V1/PID ALTURA/Sum Fdbk'
 * '<S84>'  : 'pid_control_V1/PID ALTURA/Tracking Mode'
 * '<S85>'  : 'pid_control_V1/PID ALTURA/Tracking Mode Sum'
 * '<S86>'  : 'pid_control_V1/PID ALTURA/Tsamp - Integral'
 * '<S87>'  : 'pid_control_V1/PID ALTURA/Tsamp - Ngain'
 * '<S88>'  : 'pid_control_V1/PID ALTURA/postSat Signal'
 * '<S89>'  : 'pid_control_V1/PID ALTURA/preInt Signal'
 * '<S90>'  : 'pid_control_V1/PID ALTURA/preSat Signal'
 * '<S91>'  : 'pid_control_V1/PID ALTURA/Anti-windup/Passthrough'
 * '<S92>'  : 'pid_control_V1/PID ALTURA/D Gain/Internal Parameters'
 * '<S93>'  : 'pid_control_V1/PID ALTURA/External Derivative/Error'
 * '<S94>'  : 'pid_control_V1/PID ALTURA/Filter/Cont. Filter'
 * '<S95>'  : 'pid_control_V1/PID ALTURA/Filter ICs/Internal IC - Filter'
 * '<S96>'  : 'pid_control_V1/PID ALTURA/I Gain/Internal Parameters'
 * '<S97>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain/Passthrough'
 * '<S98>'  : 'pid_control_V1/PID ALTURA/Ideal P Gain Fdbk/Disabled'
 * '<S99>'  : 'pid_control_V1/PID ALTURA/Integrator/Continuous'
 * '<S100>' : 'pid_control_V1/PID ALTURA/Integrator ICs/Internal IC'
 * '<S101>' : 'pid_control_V1/PID ALTURA/N Copy/Disabled'
 * '<S102>' : 'pid_control_V1/PID ALTURA/N Gain/Internal Parameters'
 * '<S103>' : 'pid_control_V1/PID ALTURA/P Copy/Disabled'
 * '<S104>' : 'pid_control_V1/PID ALTURA/Parallel P Gain/Internal Parameters'
 * '<S105>' : 'pid_control_V1/PID ALTURA/Reset Signal/Disabled'
 * '<S106>' : 'pid_control_V1/PID ALTURA/Saturation/Passthrough'
 * '<S107>' : 'pid_control_V1/PID ALTURA/Saturation Fdbk/Disabled'
 * '<S108>' : 'pid_control_V1/PID ALTURA/Sum/Sum_PID'
 * '<S109>' : 'pid_control_V1/PID ALTURA/Sum Fdbk/Disabled'
 * '<S110>' : 'pid_control_V1/PID ALTURA/Tracking Mode/Disabled'
 * '<S111>' : 'pid_control_V1/PID ALTURA/Tracking Mode Sum/Passthrough'
 * '<S112>' : 'pid_control_V1/PID ALTURA/Tsamp - Integral/TsSignalSpecification'
 * '<S113>' : 'pid_control_V1/PID ALTURA/Tsamp - Ngain/Passthrough'
 * '<S114>' : 'pid_control_V1/PID ALTURA/postSat Signal/Forward_Path'
 * '<S115>' : 'pid_control_V1/PID ALTURA/preInt Signal/Internal PreInt'
 * '<S116>' : 'pid_control_V1/PID ALTURA/preSat Signal/Forward_Path'
 * '<S117>' : 'pid_control_V1/PID PITCH//ELEVATOR/Anti-windup'
 * '<S118>' : 'pid_control_V1/PID PITCH//ELEVATOR/D Gain'
 * '<S119>' : 'pid_control_V1/PID PITCH//ELEVATOR/External Derivative'
 * '<S120>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter'
 * '<S121>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter ICs'
 * '<S122>' : 'pid_control_V1/PID PITCH//ELEVATOR/I Gain'
 * '<S123>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain'
 * '<S124>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain Fdbk'
 * '<S125>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator'
 * '<S126>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator ICs'
 * '<S127>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Copy'
 * '<S128>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Gain'
 * '<S129>' : 'pid_control_V1/PID PITCH//ELEVATOR/P Copy'
 * '<S130>' : 'pid_control_V1/PID PITCH//ELEVATOR/Parallel P Gain'
 * '<S131>' : 'pid_control_V1/PID PITCH//ELEVATOR/Reset Signal'
 * '<S132>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation'
 * '<S133>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation Fdbk'
 * '<S134>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum'
 * '<S135>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum Fdbk'
 * '<S136>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode'
 * '<S137>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode Sum'
 * '<S138>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Integral'
 * '<S139>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Ngain'
 * '<S140>' : 'pid_control_V1/PID PITCH//ELEVATOR/postSat Signal'
 * '<S141>' : 'pid_control_V1/PID PITCH//ELEVATOR/preInt Signal'
 * '<S142>' : 'pid_control_V1/PID PITCH//ELEVATOR/preSat Signal'
 * '<S143>' : 'pid_control_V1/PID PITCH//ELEVATOR/Anti-windup/Back Calculation'
 * '<S144>' : 'pid_control_V1/PID PITCH//ELEVATOR/D Gain/Internal Parameters'
 * '<S145>' : 'pid_control_V1/PID PITCH//ELEVATOR/External Derivative/Error'
 * '<S146>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter/Cont. Filter'
 * '<S147>' : 'pid_control_V1/PID PITCH//ELEVATOR/Filter ICs/Internal IC - Filter'
 * '<S148>' : 'pid_control_V1/PID PITCH//ELEVATOR/I Gain/Internal Parameters'
 * '<S149>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain/Passthrough'
 * '<S150>' : 'pid_control_V1/PID PITCH//ELEVATOR/Ideal P Gain Fdbk/Disabled'
 * '<S151>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator/Continuous'
 * '<S152>' : 'pid_control_V1/PID PITCH//ELEVATOR/Integrator ICs/Internal IC'
 * '<S153>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Copy/Disabled'
 * '<S154>' : 'pid_control_V1/PID PITCH//ELEVATOR/N Gain/Internal Parameters'
 * '<S155>' : 'pid_control_V1/PID PITCH//ELEVATOR/P Copy/Disabled'
 * '<S156>' : 'pid_control_V1/PID PITCH//ELEVATOR/Parallel P Gain/Internal Parameters'
 * '<S157>' : 'pid_control_V1/PID PITCH//ELEVATOR/Reset Signal/Disabled'
 * '<S158>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation/Enabled'
 * '<S159>' : 'pid_control_V1/PID PITCH//ELEVATOR/Saturation Fdbk/Disabled'
 * '<S160>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum/Sum_PID'
 * '<S161>' : 'pid_control_V1/PID PITCH//ELEVATOR/Sum Fdbk/Disabled'
 * '<S162>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode/Disabled'
 * '<S163>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tracking Mode Sum/Passthrough'
 * '<S164>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Integral/TsSignalSpecification'
 * '<S165>' : 'pid_control_V1/PID PITCH//ELEVATOR/Tsamp - Ngain/Passthrough'
 * '<S166>' : 'pid_control_V1/PID PITCH//ELEVATOR/postSat Signal/Forward_Path'
 * '<S167>' : 'pid_control_V1/PID PITCH//ELEVATOR/preInt Signal/Internal PreInt'
 * '<S168>' : 'pid_control_V1/PID PITCH//ELEVATOR/preSat Signal/Forward_Path'
 * '<S169>' : 'pid_control_V1/PID TIIMON/Anti-windup'
 * '<S170>' : 'pid_control_V1/PID TIIMON/D Gain'
 * '<S171>' : 'pid_control_V1/PID TIIMON/External Derivative'
 * '<S172>' : 'pid_control_V1/PID TIIMON/Filter'
 * '<S173>' : 'pid_control_V1/PID TIIMON/Filter ICs'
 * '<S174>' : 'pid_control_V1/PID TIIMON/I Gain'
 * '<S175>' : 'pid_control_V1/PID TIIMON/Ideal P Gain'
 * '<S176>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk'
 * '<S177>' : 'pid_control_V1/PID TIIMON/Integrator'
 * '<S178>' : 'pid_control_V1/PID TIIMON/Integrator ICs'
 * '<S179>' : 'pid_control_V1/PID TIIMON/N Copy'
 * '<S180>' : 'pid_control_V1/PID TIIMON/N Gain'
 * '<S181>' : 'pid_control_V1/PID TIIMON/P Copy'
 * '<S182>' : 'pid_control_V1/PID TIIMON/Parallel P Gain'
 * '<S183>' : 'pid_control_V1/PID TIIMON/Reset Signal'
 * '<S184>' : 'pid_control_V1/PID TIIMON/Saturation'
 * '<S185>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk'
 * '<S186>' : 'pid_control_V1/PID TIIMON/Sum'
 * '<S187>' : 'pid_control_V1/PID TIIMON/Sum Fdbk'
 * '<S188>' : 'pid_control_V1/PID TIIMON/Tracking Mode'
 * '<S189>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum'
 * '<S190>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral'
 * '<S191>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain'
 * '<S192>' : 'pid_control_V1/PID TIIMON/postSat Signal'
 * '<S193>' : 'pid_control_V1/PID TIIMON/preInt Signal'
 * '<S194>' : 'pid_control_V1/PID TIIMON/preSat Signal'
 * '<S195>' : 'pid_control_V1/PID TIIMON/Anti-windup/Passthrough'
 * '<S196>' : 'pid_control_V1/PID TIIMON/D Gain/Internal Parameters'
 * '<S197>' : 'pid_control_V1/PID TIIMON/External Derivative/Error'
 * '<S198>' : 'pid_control_V1/PID TIIMON/Filter/Cont. Filter'
 * '<S199>' : 'pid_control_V1/PID TIIMON/Filter ICs/Internal IC - Filter'
 * '<S200>' : 'pid_control_V1/PID TIIMON/I Gain/Internal Parameters'
 * '<S201>' : 'pid_control_V1/PID TIIMON/Ideal P Gain/Passthrough'
 * '<S202>' : 'pid_control_V1/PID TIIMON/Ideal P Gain Fdbk/Disabled'
 * '<S203>' : 'pid_control_V1/PID TIIMON/Integrator/Continuous'
 * '<S204>' : 'pid_control_V1/PID TIIMON/Integrator ICs/Internal IC'
 * '<S205>' : 'pid_control_V1/PID TIIMON/N Copy/Disabled'
 * '<S206>' : 'pid_control_V1/PID TIIMON/N Gain/Internal Parameters'
 * '<S207>' : 'pid_control_V1/PID TIIMON/P Copy/Disabled'
 * '<S208>' : 'pid_control_V1/PID TIIMON/Parallel P Gain/Internal Parameters'
 * '<S209>' : 'pid_control_V1/PID TIIMON/Reset Signal/Disabled'
 * '<S210>' : 'pid_control_V1/PID TIIMON/Saturation/Enabled'
 * '<S211>' : 'pid_control_V1/PID TIIMON/Saturation Fdbk/Disabled'
 * '<S212>' : 'pid_control_V1/PID TIIMON/Sum/Sum_PID'
 * '<S213>' : 'pid_control_V1/PID TIIMON/Sum Fdbk/Disabled'
 * '<S214>' : 'pid_control_V1/PID TIIMON/Tracking Mode/Disabled'
 * '<S215>' : 'pid_control_V1/PID TIIMON/Tracking Mode Sum/Passthrough'
 * '<S216>' : 'pid_control_V1/PID TIIMON/Tsamp - Integral/TsSignalSpecification'
 * '<S217>' : 'pid_control_V1/PID TIIMON/Tsamp - Ngain/Passthrough'
 * '<S218>' : 'pid_control_V1/PID TIIMON/postSat Signal/Forward_Path'
 * '<S219>' : 'pid_control_V1/PID TIIMON/preInt Signal/Internal PreInt'
 * '<S220>' : 'pid_control_V1/PID TIIMON/preSat Signal/Forward_Path'
 * '<S221>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup'
 * '<S222>' : 'pid_control_V1/PID VELOCIDAD/D Gain'
 * '<S223>' : 'pid_control_V1/PID VELOCIDAD/External Derivative'
 * '<S224>' : 'pid_control_V1/PID VELOCIDAD/Filter'
 * '<S225>' : 'pid_control_V1/PID VELOCIDAD/Filter ICs'
 * '<S226>' : 'pid_control_V1/PID VELOCIDAD/I Gain'
 * '<S227>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain'
 * '<S228>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain Fdbk'
 * '<S229>' : 'pid_control_V1/PID VELOCIDAD/Integrator'
 * '<S230>' : 'pid_control_V1/PID VELOCIDAD/Integrator ICs'
 * '<S231>' : 'pid_control_V1/PID VELOCIDAD/N Copy'
 * '<S232>' : 'pid_control_V1/PID VELOCIDAD/N Gain'
 * '<S233>' : 'pid_control_V1/PID VELOCIDAD/P Copy'
 * '<S234>' : 'pid_control_V1/PID VELOCIDAD/Parallel P Gain'
 * '<S235>' : 'pid_control_V1/PID VELOCIDAD/Reset Signal'
 * '<S236>' : 'pid_control_V1/PID VELOCIDAD/Saturation'
 * '<S237>' : 'pid_control_V1/PID VELOCIDAD/Saturation Fdbk'
 * '<S238>' : 'pid_control_V1/PID VELOCIDAD/Sum'
 * '<S239>' : 'pid_control_V1/PID VELOCIDAD/Sum Fdbk'
 * '<S240>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode'
 * '<S241>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode Sum'
 * '<S242>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Integral'
 * '<S243>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Ngain'
 * '<S244>' : 'pid_control_V1/PID VELOCIDAD/postSat Signal'
 * '<S245>' : 'pid_control_V1/PID VELOCIDAD/preInt Signal'
 * '<S246>' : 'pid_control_V1/PID VELOCIDAD/preSat Signal'
 * '<S247>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel'
 * '<S248>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel/Dead Zone'
 * '<S249>' : 'pid_control_V1/PID VELOCIDAD/Anti-windup/Cont. Clamping Parallel/Dead Zone/Enabled'
 * '<S250>' : 'pid_control_V1/PID VELOCIDAD/D Gain/Internal Parameters'
 * '<S251>' : 'pid_control_V1/PID VELOCIDAD/External Derivative/Error'
 * '<S252>' : 'pid_control_V1/PID VELOCIDAD/Filter/Cont. Filter'
 * '<S253>' : 'pid_control_V1/PID VELOCIDAD/Filter ICs/Internal IC - Filter'
 * '<S254>' : 'pid_control_V1/PID VELOCIDAD/I Gain/Internal Parameters'
 * '<S255>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain/Passthrough'
 * '<S256>' : 'pid_control_V1/PID VELOCIDAD/Ideal P Gain Fdbk/Disabled'
 * '<S257>' : 'pid_control_V1/PID VELOCIDAD/Integrator/Continuous'
 * '<S258>' : 'pid_control_V1/PID VELOCIDAD/Integrator ICs/Internal IC'
 * '<S259>' : 'pid_control_V1/PID VELOCIDAD/N Copy/Disabled'
 * '<S260>' : 'pid_control_V1/PID VELOCIDAD/N Gain/Internal Parameters'
 * '<S261>' : 'pid_control_V1/PID VELOCIDAD/P Copy/Disabled'
 * '<S262>' : 'pid_control_V1/PID VELOCIDAD/Parallel P Gain/Internal Parameters'
 * '<S263>' : 'pid_control_V1/PID VELOCIDAD/Reset Signal/Disabled'
 * '<S264>' : 'pid_control_V1/PID VELOCIDAD/Saturation/Enabled'
 * '<S265>' : 'pid_control_V1/PID VELOCIDAD/Saturation Fdbk/Disabled'
 * '<S266>' : 'pid_control_V1/PID VELOCIDAD/Sum/Sum_PID'
 * '<S267>' : 'pid_control_V1/PID VELOCIDAD/Sum Fdbk/Disabled'
 * '<S268>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode/Disabled'
 * '<S269>' : 'pid_control_V1/PID VELOCIDAD/Tracking Mode Sum/Passthrough'
 * '<S270>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Integral/TsSignalSpecification'
 * '<S271>' : 'pid_control_V1/PID VELOCIDAD/Tsamp - Ngain/Passthrough'
 * '<S272>' : 'pid_control_V1/PID VELOCIDAD/postSat Signal/Forward_Path'
 * '<S273>' : 'pid_control_V1/PID VELOCIDAD/preInt Signal/Internal PreInt'
 * '<S274>' : 'pid_control_V1/PID VELOCIDAD/preSat Signal/Forward_Path'
 * '<S275>' : 'pid_control_V1/Subscribe-YAW/Enabled Subsystem'
 * '<S276>' : 'pid_control_V1/Subscribe-YAW1/Enabled Subsystem'
 * '<S277>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))'
 * '<S278>' : 'pid_control_V1/Subsystem/MATLAB Function'
 * '<S279>' : 'pid_control_V1/Subsystem/Subscribe'
 * '<S280>' : 'pid_control_V1/Subsystem/Subscribe1'
 * '<S281>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Angle Conversion'
 * '<S282>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates'
 * '<S283>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities'
 * '<S284>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Length Conversion'
 * '<S285>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Length Conversion1'
 * '<S286>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities'
 * '<S287>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates'
 * '<S288>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities'
 * '<S289>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths'
 * '<S290>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Velocity Conversion'
 * '<S291>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Velocity Conversion2'
 * '<S292>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/White Noise'
 * '<S293>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hpgw'
 * '<S294>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hqgw'
 * '<S295>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on angular rates/Hrgw'
 * '<S296>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hugw(s)'
 * '<S297>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hvgw(s)'
 * '<S298>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Filters on velocities/Hwgw(s)'
 * '<S299>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities/High Altitude Intensity'
 * '<S300>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/RMS turbulence  intensities/Low Altitude Intensity'
 * '<S301>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates'
 * '<S302>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates'
 * '<S303>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Medium//High  altitude rates'
 * '<S304>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Merge Subsystems'
 * '<S305>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates/wind to body transformation'
 * '<S306>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Interpolate  rates/wind to body transformation/convert to earth coords'
 * '<S307>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates/wind to body transformation'
 * '<S308>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select angular rates/Low altitude  rates/wind to body transformation/convert to earth coords'
 * '<S309>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities'
 * '<S310>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities'
 * '<S311>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Medium//High  altitude velocities'
 * '<S312>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Merge Subsystems'
 * '<S313>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities/wind to body transformation'
 * '<S314>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Interpolate  velocities/wind to body transformation/convert to earth coords'
 * '<S315>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities/wind to body transformation'
 * '<S316>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Select velocities/Low altitude  velocities/wind to body transformation/convert to earth coords'
 * '<S317>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Low altitude scale length'
 * '<S318>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Medium//High altitude scale length'
 * '<S319>' : 'pid_control_V1/Subsystem/Dryden Wind Turbulence Model  (Continuous (-q +r))/Turbulence scale lengths/Medium//High altitude scale length/Length Conversion'
 * '<S320>' : 'pid_control_V1/Subsystem/Subscribe/Enabled Subsystem'
 * '<S321>' : 'pid_control_V1/Subsystem/Subscribe1/Enabled Subsystem'
 */
#endif                                 /* pid_control_V1_h_ */
