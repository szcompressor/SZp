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
#include "szp_utility.h"
#include "szp_defines.h"
#include "szp_float.h"
#include "szpd_float.h"
#include "szp_TypeManager.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "szp_defines.h"
#include "szp_float.h"
#include "szpd_float.h"
#include "szp_TypeManager.h"

typedef union lint16
{
	unsigned short usvalue;
	short svalue;
	unsigned char byte[2];
} lint16;

typedef union lint32
{
	int ivalue;
	unsigned int uivalue;
	unsigned char byte[4];
} lint32;

typedef union lint64
{
	long lvalue;
	unsigned long ulvalue;
	unsigned char byte[8];
} lint64;

typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;


extern int versionNumber[4];

//-------------------key global variables--------------
extern int dataEndianType; //*endian type of the data read from disk
extern int sysEndianType; //*sysEndianType is actually set automatically.

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szp_H  ----- */
