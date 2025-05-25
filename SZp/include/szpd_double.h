/**
 *  @file szpd_double.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szpd_double_H
#define _szpd_double_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include "szp_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

void szp_double_decompress_openmp_threadblock(double **newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_double_decompress_openmp_threadblock_arg(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_double_decompress_single_thread_arg(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes);

size_t szp_double_decompress_single_thread_arg_record(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes);

void szp_double_decompress_openmp_threadblock_randomaccess_args(double **newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szpd_double_H  ----- */
