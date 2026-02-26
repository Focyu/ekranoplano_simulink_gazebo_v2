#ifndef OPEN_LOOP_V1__VISIBILITY_CONTROL_H_
#define OPEN_LOOP_V1__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define OPEN_LOOP_V1_EXPORT __attribute__ ((dllexport))
    #define OPEN_LOOP_V1_IMPORT __attribute__ ((dllimport))
  #else
    #define OPEN_LOOP_V1_EXPORT __declspec(dllexport)
    #define OPEN_LOOP_V1_IMPORT __declspec(dllimport)
  #endif
  #ifdef OPEN_LOOP_V1_BUILDING_LIBRARY
    #define OPEN_LOOP_V1_PUBLIC OPEN_LOOP_V1_EXPORT
  #else
    #define OPEN_LOOP_V1_PUBLIC OPEN_LOOP_V1_IMPORT
  #endif
  #define OPEN_LOOP_V1_PUBLIC_TYPE OPEN_LOOP_V1_PUBLIC
  #define OPEN_LOOP_V1_LOCAL
#else
  #define OPEN_LOOP_V1_EXPORT __attribute__ ((visibility("default")))
  #define OPEN_LOOP_V1_IMPORT
  #if __GNUC__ >= 4
    #define OPEN_LOOP_V1_PUBLIC __attribute__ ((visibility("default")))
    #define OPEN_LOOP_V1_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define OPEN_LOOP_V1_PUBLIC
    #define OPEN_LOOP_V1_LOCAL
  #endif
  #define OPEN_LOOP_V1_PUBLIC_TYPE
#endif
#endif  // OPEN_LOOP_V1__VISIBILITY_CONTROL_H_
// Generated 26-Feb-2026 18:13:34
 