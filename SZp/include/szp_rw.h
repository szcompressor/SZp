/**
 *  @file szp_rw.h
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief Header file for the whole io interface.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _szp_RW_H
#define _szp_RW_H

#include <stdio.h>
#include <stdint.h>
#include "szp_CompressionToolkit.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

int szp_checkFileExistance(char* filePath);

float** szp_create2DArray_float(size_t m, size_t n);
void szp_free2DArray_float(float** data, size_t m);
float*** szp_create3DArray_float(size_t p, size_t m, size_t n);
void szp_free3DArray_float(float*** data, size_t p, size_t m);
double** szp_create2DArray_double(size_t m, size_t n);
void szp_free2DArray_double(double** data, size_t m);
double*** szp_create3DArray_double(size_t p, size_t m, size_t n);
void szp_free3DArray_double(double*** data, size_t p, size_t m);
size_t szp_checkFileSize(char *srcFilePath, int *status);

unsigned char *szp_readByteData(char *srcFilePath, size_t *byteLength, int *status);
double *szp_readDoubleData(char *srcFilePath, size_t *nbEle, int *status);
int8_t *szp_readInt8Data(char *srcFilePath, size_t *nbEle, int *status);
int16_t *szp_readInt16Data(char *srcFilePath, size_t *nbEle, int *status);
uint16_t *szp_readUInt16Data(char *srcFilePath, size_t *nbEle, int *status);
int32_t *szp_readInt32Data(char *srcFilePath, size_t *nbEle, int *status);
uint32_t *szp_readUInt32Data(char *srcFilePath, size_t *nbEle, int *status);
int64_t *szp_readInt64Data(char *srcFilePath, size_t *nbEle, int *status);
uint64_t *szp_readUInt64Data(char *srcFilePath, size_t *nbEle, int *status);
float *szp_readFloatData(char *srcFilePath, size_t *nbEle, int *status);
unsigned short* szp_readShortData(char *srcFilePath, size_t *dataLength, int *status);

double *szp_readDoubleData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
int8_t *szp_readInt8Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
int16_t *szp_readInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
uint16_t *szp_readUInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
int32_t *szp_readInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
uint32_t *szp_readUInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
int64_t *szp_readInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
uint64_t *szp_readUInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
float *szp_readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);

void szp_writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status);
void szp_writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status);
void szp_writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status);
void szp_writeData(void *data, int dataType, size_t nbEle, char *tgtFilePath, int *status);
void szp_writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status);
void szp_writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status);
void szp_writeShortData_inBytes(short *states, size_t stateLength, char *tgtFilePath, int *status);
void szp_writeUShortData_inBytes(unsigned short *states, size_t stateLength, char *tgtFilePath, int *status);
void szp_writeIntData_inBytes(int *states, size_t stateLength, char *tgtFilePath, int *status);
void szp_writeUIntData_inBytes(unsigned int *states, size_t stateLength, char *tgtFilePath, int *status);
void szp_writeLongData_inBytes(int64_t *states, size_t stateLength, char *tgtFilePath, int *status);
void szp_writeULongData_inBytes(uint64_t *states, size_t stateLength, char *tgtFilePath, int *status);

void szp_writeStrings(int nbStr, char *str[], char *tgtFilePath, int *status);

//void szp_convertToPFM_float(float *data, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, int endianType, char *tgtFilePath, int *status);

void szp_checkfilesizec_(char *srcFilePath, int *len, size_t *filesize);
void szp_readbytefile_(char *srcFilePath, int *len, unsigned char *bytes, size_t *byteLength);
void szp_readdoublefile_(char *srcFilePath, int *len, double *data, size_t *nbEle);
void szp_readfloatfile_(char *srcFilePath, int *len, float *data, size_t *nbEle);
void szp_writebytefile_(unsigned char *bytes, size_t *byteLength, char *tgtFilePath, int *len);
void szp_writedoublefile_(double *data, size_t *nbEle, char *tgtFilePath, int *len);
void szp_writefloatfile_(float *data, size_t *nbEle, char *tgtFilePath, int *len);

}

#endif /* ----- #ifndef _szp_RW_H  ----- */
