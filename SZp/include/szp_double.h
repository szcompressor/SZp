/**
 *  @file szp_double.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szp_double_H
#define _szp_double_H

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
szp_double_openmp_direct_predict_quantization(double *oriData, size_t *outSize, double absErrBound,
                                             size_t nbEle, int blockSize);

int *
szp_double_openmp_threadblock_predict_quantization(double *oriData, size_t *outSize, double absErrBound,
                                                  size_t nbEle, int blockSize);

unsigned char *
szp_double_openmp_threadblock(double *oriData, size_t *outSize, double absErrBound,
                             size_t nbEle, int blockSize);

void szp_double_openmp_threadblock_arg(unsigned char *output, double *oriData, size_t *outSize, double absErrBound,
                                      size_t nbEle, int blockSize);

void szp_double_single_thread_arg(unsigned char *output, double *oriData, size_t *outSize, double absErrBound,
                                 size_t nbEle, int blockSize);

size_t szp_double_single_thread_arg_record(unsigned char *output, double *oriData, size_t *outSize, double absErrBound,
                                       size_t nbEle, int blockSize);

void
szp_double_openmp_threadblock_randomaccess_args(unsigned char* output, double *oriData, size_t *outSize, double absErrBound,
                                          size_t nbEle, int blockSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szp_double_H  ----- */
