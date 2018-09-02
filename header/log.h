#ifndef _LOG_H_
#define _LOG_H_

//
//  log.h
//
//  Created by Jaewook Kang
//

/* ----------- Define OS TYPE or Simulator ---------- */
// __ANDROID__				// (1) Android OS
// __APPLE__				// (2) iOS
/* -------------------------------------------------- */

/* ----------------- (1) Android OS ----------------- */
#ifdef __ANDROID__

// Data logging
#undef MODEM_DATA_LOG
#undef PACKET_COMBINER_LOG

#ifndef __JENKINSBUILD__
#include <android/log.h>
#endif

#define  LOG_TAG "EveryDSP"

// Common Logs for All BuildTypes : Error / Warn
#define  LOGE(...) __android_log_print(ANDROID_LOG_ERROR  ,LOG_TAG,__VA_ARGS__)
#define  LOGW(...) __android_log_print(ANDROID_LOG_WARN   ,LOG_TAG,__VA_ARGS__)

// Conditional Logs : Info / Debug / Verbose
#ifdef _LOG_DEBUG
#define  LOGI(...) __android_log_print(ANDROID_LOG_INFO   ,LOG_TAG,__VA_ARGS__)
#define  LOGD(...) __android_log_print(ANDROID_LOG_DEBUG  ,LOG_TAG,__VA_ARGS__)
#define  LOGV(...) __android_log_print(ANDROID_LOG_VERBOSE,LOG_TAG,__VA_ARGS__)
#else
#define  LOGI(...)
#define  LOGD(...)
#define  LOGV(...)
#endif

#endif
/* -------------------------------------------------- */

/* --------------------- (2) iOS -------------------- */
// iOS LOG macro ..
#ifdef __APPLE__
#include <stdio.h>
#include <iostream>
//void logi(std::string format, ...);
//void loge(std::string format, ...);
//void logv(std::string format, ...);
//void logw(std::string format, ...);
//void logd(std::string format, ...);
//#define  LOGE(...) loge(__VA_ARGS__)
//#define  LOGW(...) logw(__VA_ARGS__)
////#define  LOGI(...) logi(__VA_ARGS__)
////#define  LOGD(...) logd(__VA_ARGS__)
//#define  LOGI(...) logv(__VA_ARGS__)
//#define  LOGD(...) logv(__VA_ARGS__)
//#define  LOGV(...) logv(__VA_ARGS__)
//
#ifdef DEBUG
#define  LOGE(...) {printf("MEX_LOG_ERROR:       "); printf(__VA_ARGS__); printf("\n");}
#define  LOGW(...) {printf("MEX_LOG_WARN:        "); printf(__VA_ARGS__); printf("\n");}
#define  LOGI(...) {printf("MEX_LOG_INFO:        "); printf(__VA_ARGS__); printf("\n");}
#define  LOGD(...) {printf("MEX_LOG_DEBUG:       "); printf(__VA_ARGS__); printf("\n");}
#define  LOGV(...) {printf("MEX_LOG_VERBOSE:     "); printf(__VA_ARGS__); printf("\n");}
#else
#define  LOGE(...) {printf("MEX_LOG_ERROR:       "); printf(__VA_ARGS__); printf("\n");}
#define  LOGW(...) {printf("MEX_LOG_WARN:        "); printf(__VA_ARGS__); printf("\n");}
#define  LOGI(...) {printf("MEX_LOG_INFO:        "); printf(__VA_ARGS__); printf("\n");}
#define  LOGD(...)
#define  LOGV(...)
#endif

#endif
/* -------------------------------------------------- */

#endif // _LOG_H_
