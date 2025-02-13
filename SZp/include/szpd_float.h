/**
 *  @file szpd_Float.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szpd_Float_H
#define _szpd_Float_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include "szp_defines.h"

void szp_float_decompress_openmp_threadblock(float **newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_float_decompress_openmp_threadblock_arg(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_float_decompress_single_thread_arg(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes);

size_t szp_float_decompress_single_thread_arg_record(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_float_decompress_openmp_threadblock_randomaccess(float **newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes);

#endif /* ----- #ifndef _szpd_Float_H  ----- */