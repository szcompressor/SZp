/**
 *  @file szp_defines.h
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief Header file for the dataCompression.c.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _szp_DEFINES_H
#define _szp_DEFINES_H

#define szp_VERNUM 0x0200
#define szp_VER_MAJOR 1
#define szp_VER_MINOR 0
#define szp_VER_BUILD 0
#define szp_VER_REVISION 0

#define ABS 0
#define REL 1
#define VR_REL 1  //alternative name to REL
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4
#define NORM 5

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14

#define SZP_NONRANDOMACCESS 1
#define SZP_RANDOMACCESS 2

#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9

#define LITTLE_ENDIAN_DATA 0 //refers to the endian type of the data read from the disk
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0 //refers to the endian type of the system
#define BIG_ENDIAN_SYSTEM 1

#define szp_NO_BLOCK_FAST_CMPR 1
#define szp_WITH_BLOCK_FAST_CMPR 2
#define szp_RANDOMACCESS_FAST_CMPR 3
#define szp_OPENMP_FAST_CMPR 4

//SUCCESS returning status
#define SZ_SCES 0  //successful
#define SZ_NSCS -1 //Not successful
#define SZ_FERR -2 //Failed to open input file
#define SZ_TERR -3 //wrong data type (should be only float or double)
#define SZ_DERR -4 //dimension error
#define SZ_MERR -5 //sz_mode error
#define SZ_BERR -6 //bound-mode error (should be only ABS, REL, ABS_AND_REL, ABS_OR_REL, or PW_REL)

typedef union szp_lint16
{
	unsigned short usvalue;
	short svalue;
	unsigned char byte[2];
} szp_lint16;

typedef union szp_lint32
{
	int ivalue;
	unsigned int uivalue;
	unsigned char byte[4];
} szp_lint32;

typedef union szp_lint64
{
	long lvalue;
	unsigned long ulvalue;
	unsigned char byte[8];
} szp_lint64;

typedef union szp_ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} szp_ldouble;

typedef union szp_lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} szp_lfloat;


extern int szp_versionNumber[4];

//-------------------key global variables--------------
extern int szp_dataEndianType; //*endian type of the data read from disk
extern int szp_sysEndianType; //*sysEndianType is actually set automatically.

#ifdef __cplusplus
extern "C" {
#endif

inline size_t szp_computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0)
	{
		dataLength = 0;
	}
	else if(r2==0)
	{
		dataLength = r1;
	}
	else if(r3==0)
	{
		dataLength = r1*r2;
	}
	else if(r4==0)
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0)
	{
		dataLength = r1*r2*r3*r4;
	}
	else
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

#ifdef __cplusplus
}
#endif

#endif /* _szp_DEFINES_H */
