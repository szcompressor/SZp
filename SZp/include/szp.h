/**
 *  @file szp.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szp_H
#define _szp_H

#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */
#include <math.h>
#include "szp_rw.h"
#include "szp_defines.h"
#include "szp_float.h"
#include "szpd_float.h"
#include "szp_double.h"
#include "szpd_double.h"
#include "szp_TypeManager.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

unsigned char *szp_fast_compress_args(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, float absErrBound,
                                      float relBoundRatio, size_t nbEle, int blockSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szp_H  ----- */
