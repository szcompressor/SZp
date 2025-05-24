/**
 *  @file szp.c
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief 
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "szp.h"
#include "szp_rw.h"

int versionNumber[4] = {szp_VER_MAJOR,szp_VER_MINOR,szp_VER_BUILD,szp_VER_REVISION};

int dataEndianType = LITTLE_ENDIAN_DATA; //*endian type of the data read from disk
int sysEndianType = LITTLE_ENDIAN_SYSTEM; //*sysEndianType is actually set automatically.

using namespace szp;

unsigned char *szp_fast_compress_args(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, float absErrBound,
                                      float relBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    unsigned char *bytes = NULL;
    size_t length = szp_computeDataLength(r5, r4, r3, r2, r1);
    size_t i = 0;
    int blockSize = 128;
    if (dataType == SZ_FLOAT)
    {

        float realPrecision = absErrBound;
        if (errBoundMode == REL)
        {
            float *oriData = (float *)data;
            float min = oriData[0];
            float max = oriData[0];
            for (i = 0; i < length; i++)
            {
                float v = oriData[i];
                if (min > v)
                    min = v;
                else if (max < v)
                    max = v;
            }
            float valueRange = max - min;
            realPrecision = valueRange * relBoundRatio;
            //printf("REAL ERROR BOUND IS %20f\n", realPrecision);
            if (fastMode == 1)
            {
                bytes = szp_float_openmp_threadblock(oriData, outSize, realPrecision,
                                                     length, blockSize);
            }
            else if (fastMode == 2)
            {
                bytes = szp_float_openmp_threadblock_randomaccess(oriData, outSize, realPrecision,
                                                                  length, blockSize);
            }
        }
        if (fastMode == 1)
        {
            bytes = szp_float_openmp_threadblock((float *)data, outSize, realPrecision,
                                                 length, blockSize);
        }
        else if (fastMode == 2)
        {
            bytes = szp_float_openmp_threadblock_randomaccess((float *)data, outSize, realPrecision,
                                                              length, blockSize);
        }
    }
    else
    {
        printf("szp only supports float type for now\n");
    }
    return bytes;
}
