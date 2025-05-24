/**
 *  @file szp_Float.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szp_Float_H
#define _szp_Float_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include "szp_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

int *
szp_float_openmp_direct_predict_quantization(float *oriData, size_t *outSize, float absErrBound,
                                             size_t nbEle, int blockSize);

int *
szp_float_openmp_threadblock_predict_quantization(float *oriData, size_t *outSize, float absErrBound,
                                                  size_t nbEle, int blockSize);

unsigned char *
szp_float_openmp_threadblock(float *oriData, size_t *outSize, float absErrBound,
                             size_t nbEle, int blockSize);

void szp_float_openmp_threadblock_arg(unsigned char *outputBytes, float *oriData, size_t *outSize, float absErrBound,
                                      size_t nbEle, int blockSize);

void szp_float_single_thread_arg(unsigned char *outputBytes, float *oriData, size_t *outSize, float absErrBound,
                                 size_t nbEle, int blockSize);

size_t szp_float_single_thread_arg_record(unsigned char *outputBytes, float *oriData, size_t *outSize, float absErrBound,
                                       size_t nbEle, int blockSize);

unsigned char *
szp_float_openmp_threadblock_randomaccess(float *oriData, size_t *outSize, float absErrBound,
                                          size_t nbEle, int blockSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szp_Float_H  ----- */
