/**
 *  @file szp_TypeManager.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#ifndef _szp_TypeManager_H
#define _szp_TypeManager_H

#include <stdio.h>
#include <stdint.h>

namespace szp{

size_t Jiajun_save_fixed_length_bits(unsigned int *unsignintArray, size_t intArrayLength, unsigned char *result, unsigned int bit_count);
size_t Jiajun_convertUInt2Byte_fast_1b_args(unsigned int *intArray, size_t intArrayLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_2b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_3b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_4b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_5b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_6b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_convertUInt2Byte_fast_7b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t Jiajun_extract_fixed_length_bits(unsigned char *result, size_t intArrayLength, unsigned int *unsignintArray, unsigned int bit_count);
void Jiajun_convertByte2UInt_fast_1b_args(size_t intArrayLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_2b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_3b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_4b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_5b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_6b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);
void Jiajun_convertByte2UInt_fast_7b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray);



size_t convertIntArray2ByteArray_fast_1b_args(unsigned char* intArray, size_t intArrayLength, unsigned char *result);
size_t convertIntArray2ByteArray_fast_1b(unsigned char* intArray, size_t intArrayLength, unsigned char **result);
size_t convertIntArray2ByteArray_fast_1b_to_result(unsigned char* intArray, size_t intArrayLength, unsigned char *result);
void convertByteArray2IntArray_fast_1b_args(size_t intArrayLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char* intArray);
void convertByteArray2IntArray_fast_1b(size_t intArrayLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_2b_args(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char *result);
size_t convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_3b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_3b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_4b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_4b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_5b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_5b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_6b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_6b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
size_t convertIntArray2ByteArray_fast_7b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_7b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);
int getLeftMovingSteps(size_t k, unsigned char resiBitLength);

}

#endif /* ----- #ifndef _szp_TypeManager_H  ----- */

