/**
 *  @file szp.cc
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *    See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "szp.h"
#include "szp_rw.h"

int szp_versionNumber[4] = {szp_VER_MAJOR, szp_VER_MINOR, szp_VER_BUILD, szp_VER_REVISION};

int szp_dataEndianType = LITTLE_ENDIAN_DATA;
int szp_sysEndianType = LITTLE_ENDIAN_SYSTEM;
using namespace szp;

unsigned char *szp_compress(int fastMode, int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound,
                            double relBoundRatio, size_t nbEle, int blockSize)
{
    unsigned char *bytes = NULL;
    size_t length = nbEle;
    size_t i = 0;
    if (dataType == SZ_FLOAT)
    {
        // For float data, we cast the double precision down to float. This is the correct behavior.
        float realPrecision = (float)absErrBound;
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
            realPrecision = valueRange * (float)relBoundRatio;
        }

        float *oriData = (float *)data;
        if (fastMode == SZP_NONRANDOMACCESS)
        {
            // bytes = (unsigned char *)malloc(sizeof(float) + sizeof(float) * length);
            // szp_float_openmp_threadblock_arg(bytes, oriData, outSize, realPrecision,
            //                                  length, blockSize); 

            // --- NON _arg version ---
            bytes = szp_float_openmp_threadblock(oriData, outSize, realPrecision, length, blockSize);
            // --- NON _arg version ---
        }
        else if (fastMode == SZP_RANDOMACCESS)
        {
            // bytes = (unsigned char *)malloc(sizeof(float) + sizeof(float) * length);
            // szp_float_openmp_threadblock_randomaccess_arg(bytes, oriData, outSize, realPrecision, length, blockSize); 

            // --- NON _arg version ---
            bytes = szp_float_openmp_threadblock_randomaccess(oriData, outSize, realPrecision, length, blockSize);  
            // --- NON _arg version ---
        }
    }
    else // SZ_DOUBLE
    {
        // For double data, we now use the full double precision from the interface.
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
            // This is now a full precision double-precision multiplication
            realPrecision = valueRange * relBoundRatio;
        }

        double *oriData = (double *)data;
        if (fastMode == SZP_NONRANDOMACCESS)
        {
            // bytes = (unsigned char *)malloc(sizeof(double) + sizeof(double) * length);
            // szp_double_openmp_threadblock_arg(bytes, oriData, outSize, realPrecision,
            //                                   length, blockSize);  

            // --- NON _arg version ---
            bytes = szp_double_openmp_threadblock(oriData, outSize, realPrecision, length, blockSize); 
            // --- NON _arg version ---
        }
        else if (fastMode == SZP_RANDOMACCESS)
        {
            // bytes = (unsigned char *)malloc(sizeof(double) + sizeof(double) * length);
            // szp_double_openmp_threadblock_randomaccess_arg(bytes, oriData, outSize, realPrecision,
            //                                                length, blockSize); 

            // --- NON _arg version ---
            bytes = szp_double_openmp_threadblock_randomaccess(oriData, outSize, realPrecision, length, blockSize); 
            // --- NON _arg version ---
        }
    }
    return bytes;
}

void *szp_decompress(int fastMode, int dataType, unsigned char *bytes, size_t byteLength, size_t nbEle, int blockSize)
{
    int x = 1;
    char *y = (char *)&x;
    if (*y == 1)
        szp_sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
        szp_sysEndianType = BIG_ENDIAN_SYSTEM;

    if (dataType == SZ_FLOAT)
    {
        
        float absErrBound = bytesToFloat(bytes);
        unsigned char *cmpBytes = bytes + sizeof(float);
        float *decompressedData = NULL;
        if (fastMode == SZP_NONRANDOMACCESS)
        {
            // decompressedData = (float *)malloc(nbEle * sizeof(float));
            // szp_float_decompress_openmp_threadblock_arg(decompressedData, nbEle, absErrBound, blockSize, cmpBytes); 

            // --- NON _arg version ---
            decompressedData = szp_float_decompress_openmp_threadblock(nbEle, absErrBound, blockSize, cmpBytes); 
            // --- NON _arg version ---
        }
        else // SZP_RANDOMACCESS
        {
            // decompressedData = (float *)malloc(nbEle * sizeof(float));
            // szp_float_decompress_openmp_threadblock_randomaccess_arg(decompressedData, nbEle, absErrBound, blockSize, cmpBytes); 

            // --- NON _arg version ---
            decompressedData = szp_float_decompress_openmp_threadblock_randomaccess(nbEle, absErrBound, blockSize, cmpBytes); 
            // --- NON _arg version ---

        }
        return decompressedData;
    }
    else // SZ_DOUBLE
    {
        double *decompressedData = NULL;
        double absErrBound = bytesToDouble(bytes);
        unsigned char *cmpBytes = bytes + sizeof(double);

        if (fastMode == SZP_NONRANDOMACCESS)
        {
            // decompressedData = (double *)malloc(nbEle * sizeof(double));
            // szp_double_decompress_openmp_threadblock_arg(decompressedData, nbEle, absErrBound, blockSize, cmpBytes); 

            // --- NON _arg version ---
            decompressedData = szp_double_decompress_openmp_threadblock(nbEle, absErrBound, blockSize, cmpBytes); 
            // --- NON _arg version ---
        }
        else // SZP_RANDOMACCESS
        {
            // decompressedData = (double *)malloc(nbEle * sizeof(double));
            // szp_double_decompress_openmp_threadblock_randomaccess_arg(decompressedData, nbEle, absErrBound, blockSize, cmpBytes); 
            // --- NON _arg version ---
            decompressedData = szp_double_decompress_openmp_threadblock_randomaccess(nbEle, absErrBound, blockSize, cmpBytes); 
            // --- NON _arg version ---
        }
        return decompressedData;
    }
}
