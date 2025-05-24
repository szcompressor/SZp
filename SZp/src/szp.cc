/**
 *  @file szp.cc
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
                                      float relBoundRatio, size_t nbEle, int blockSize)
{
    unsigned char *bytes = NULL;
    size_t length = nbEle;
    size_t i = 0;
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
            if (fastMode == SZP_NONRANDOMACCESS)
            {
                bytes = szp_float_openmp_threadblock(oriData, outSize, realPrecision,
                                                     length, blockSize);
            }
            else if (fastMode == SZP_RANDOMACCESS)
            {
                bytes = szp_float_openmp_threadblock_randomaccess(oriData, outSize, realPrecision,
                                                                  length, blockSize);
            }
        }
        if (fastMode == SZP_NONRANDOMACCESS)
        {
			float *oriData = (float *)data;
            bytes = szp_float_openmp_threadblock(oriData, outSize, realPrecision,
                                                 length, blockSize);
        }
        else if (fastMode == SZP_RANDOMACCESS)
        {
			float *oriData = (float *)data;
            bytes = szp_float_openmp_threadblock_randomaccess(oriData, outSize, realPrecision,
                                                              length, blockSize);
        }
    }
    else
    {

        double realPrecision = absErrBound;
        if (errBoundMode == REL)
        {
            double *oriData = (double *)data;
            double min = oriData[0];
            double max = oriData[0];
            for (i = 0; i < length; i++)
            {
                double v = oriData[i];
                if (min > v)
                    min = v;
                else if (max < v)
                    max = v;
            }
            double valueRange = max - min;
            realPrecision = valueRange * relBoundRatio;
            //printf("REAL ERROR BOUND IS %20f\n", realPrecision);
            if (fastMode == SZP_NONRANDOMACCESS)
            {
                bytes = szp_double_openmp_threadblock(oriData, outSize, realPrecision,
                                                     length, blockSize);
            }
            else if (fastMode == SZP_RANDOMACCESS)
            {
                bytes = szp_double_openmp_threadblock_randomaccess(oriData, outSize, realPrecision,
                                                                  length, blockSize);
            }
        }
        if (fastMode == SZP_NONRANDOMACCESS)
        {
			double *oriData = (double *)data;
            bytes = szp_double_openmp_threadblock(oriData, outSize, realPrecision,
                                                 length, blockSize);
        }
        else if (fastMode == SZP_RANDOMACCESS)
        {
			double *oriData = (double *)data;
            bytes = szp_double_openmp_threadblock_randomaccess(oriData, outSize, realPrecision,
                                                              length, blockSize);
        }        
    }
    return bytes;
}
