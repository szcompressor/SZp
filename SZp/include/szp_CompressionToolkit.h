/**
 *  @file szp_CompressionToolkit.h
 *  @author Sheng Di
 *  @date Feb, 2025
 *  @brief Header file for the ByteToolkit.c.
 *  (C) 2025 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _szp_CompressionToolkit_H
#define _szp_CompressionToolkit_H

template <typename T>
void computeMinMax(const T* data, size_t len, T& min, T& max) {
    min = data[0];
    max = data[0];
    for (size_t i = 1; i < len; ++i) {
        if (data[i] < min) min = data[i];
        else if (data[i] > max) max = data[i];
    }
}


#ifdef __cplusplus
extern "C" {
#endif

#include "szp_defines.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <cstdint>
        
inline void symTransform_8bytes(unsigned char data[8])
{
	unsigned char tmp = data[0];
	data[0] = data[7];
	data[7] = tmp;

	tmp = data[1];
	data[1] = data[6];
	data[6] = tmp;
	
	tmp = data[2];
	data[2] = data[5];
	data[5] = tmp;
	
	tmp = data[3];
	data[3] = data[4];
	data[4] = tmp;
}

inline void symTransform_2bytes(unsigned char data[2])
{
	unsigned char tmp = data[0];
	data[0] = data[1];
	data[1] = tmp;
}

inline void symTransform_4bytes(unsigned char data[4])
{
	unsigned char tmp = data[0];
	data[0] = data[3];
	data[3] = tmp;

	tmp = data[1];
	data[1] = data[2];
	data[2] = tmp;
}

inline void sz_writeBits_Fast_int8(unsigned char* buffer,uint64_t *bitPosPtr, int numBits, unsigned char data)
{
    unsigned char mask = (1 << numBits)-1;
    *(buffer + ((*bitPosPtr)>>3)) |= (data & mask) << ((*bitPosPtr) & (uint64_t)0x0000000000000007);
    (*bitPosPtr) += numBits;
}

inline void sz_writeBits_Fast_int32(unsigned char* buffer,uint64_t *bitPosPtr, int numBits, int32_t data)
{
    uint32_t mask = (1 << numBits)-1;
    *(uint32_t*)(buffer + ((*bitPosPtr)>>3)) |= ((*(uint32_t*)&data)&mask) << ((*bitPosPtr) & (uint64_t)0x0000000000000007);
    (*bitPosPtr) += numBits;
}

inline void sz_writeBits_Fast_int64(unsigned char* buffer,uint64_t *bitPosPtr, int numBits, int64_t data)
{
    uint64_t mask = ((uint64_t)0x0000000000000001<<numBits)-1;
    *(uint64_t*)(buffer + ((*bitPosPtr)>>3)) |= ((*(uint64_t*)&data)&mask) << ((*bitPosPtr) & (uint64_t)0x0000000000000007);
    (*bitPosPtr) += numBits;
}


inline unsigned short bytesToUInt16_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	unsigned short res = 0;
	
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;
	
	return res;
}	
	
inline unsigned int bytesToUInt32_bigEndian(unsigned char* bytes)
{
	unsigned int temp = 0;
	unsigned int res = 0;
	
	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;
	
	return res;
}

inline unsigned long bytesToUInt64_bigEndian(unsigned char* b) {
	unsigned long temp = 0;
	unsigned long res = 0;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}
	
inline short bytesToInt16_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	short res = 0;
	
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;
	
	return res;
}	
	
inline int bytesToInt32_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	int res = 0;
	
	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;
	
	return res;
}

inline long bytesToInt64_bigEndian(unsigned char* b) {
	long temp = 0;
	long res = 0;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}

inline int bytesToInt_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	int res = 0;
	
	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;
	
	return res;
}

/**
 * @unsigned char *b the variable to store the converted bytes (length=4)
 * @unsigned int num
 * */
inline void intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
	b[0] = (unsigned char)(num >> 24);	
	b[1] = (unsigned char)(num >> 16);	
	b[2] = (unsigned char)(num >> 8);	
	b[3] = (unsigned char)(num);	
	
	//note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

inline void int64ToBytes_bigEndian(unsigned char *b, uint64_t num)
{
	b[0] = (unsigned char)(num>>56);
	b[1] = (unsigned char)(num>>48);
	b[2] = (unsigned char)(num>>40);
	b[3] = (unsigned char)(num>>32);
	b[4] = (unsigned char)(num>>24);
	b[5] = (unsigned char)(num>>16);
	b[6] = (unsigned char)(num>>8);
	b[7] = (unsigned char)(num);
}

inline void int32ToBytes_bigEndian(unsigned char *b, uint32_t num)
{
	b[0] = (unsigned char)(num >> 24);	
	b[1] = (unsigned char)(num >> 16);	
	b[2] = (unsigned char)(num >> 8);	
	b[3] = (unsigned char)(num);		
}

inline void int16ToBytes_bigEndian(unsigned char *b, uint16_t num)
{
	b[0] = (unsigned char)(num >> 8);	
	b[1] = (unsigned char)(num);
}

/**
 * @endianType: refers to the endian_type of unsigned char* b.
 * */
inline long bytesToLong_bigEndian(unsigned char* b) {
	long temp = 0;
	long res = 0;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}

inline void longToBytes_bigEndian(unsigned char *b, unsigned long num) 
{
	b[0] = (unsigned char)(num>>56);
	b[1] = (unsigned char)(num>>48);
	b[2] = (unsigned char)(num>>40);
	b[3] = (unsigned char)(num>>32);
	b[4] = (unsigned char)(num>>24);
	b[5] = (unsigned char)(num>>16);
	b[6] = (unsigned char)(num>>8);
	b[7] = (unsigned char)(num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		symTransform_8bytes(*b);
}


inline long doubleToOSEndianLong(double value)
{
	ldouble buf;
	buf.value = value;
	return buf.lvalue;
}

inline int floatToOSEndianInt(float value)
{
	lfloat buf;
	buf.value = value;
	return buf.ivalue;
}

//TODO: debug: lfBuf.lvalue could be actually little_endian....
inline short getExponent_float(float value)
{
	//int ivalue = floatToBigEndianInt(value);

	lfloat lbuf;
	lbuf.value = value;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
	return (short)expValue;
}

inline short getPrecisionReqLength_float(float precision)
{
	lfloat lbuf;
	lbuf.value = precision;
	int ivalue = lbuf.ivalue;
	
	int expValue = (ivalue & 0x7F800000) >> 23;
	expValue -= 127;
//	unsigned char the1stManBit = (unsigned char)((ivalue & 0x00400000) >> 22);
//	if(the1stManBit==1)
//		expValue--;	
	return (short)expValue;
}

inline short getExponent_double(double value)
{
	//long lvalue = doubleToBigEndianLong(value);
	
	ldouble lbuf;
	lbuf.value = value;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
	return (short)expValue;
}

inline short getPrecisionReqLength_double(double precision)
{
	ldouble lbuf;
	lbuf.value = precision;
	long lvalue = lbuf.lvalue;
	
	int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
	expValue -= 1023;
//	unsigned char the1stManBit = (unsigned char)((lvalue & 0x0008000000000000) >> 51);
//	if(the1stManBit==1)
//		expValue--;
	return (short)expValue;
}

inline unsigned char numberOfLeadingZeros_Int(int i) {
	if (i == 0)
		return 32;
	unsigned char n = 1;
	if (((unsigned int)i) >> 16 == 0) { n += 16; i <<= 16; }
	if (((unsigned int)i) >> 24 == 0) { n +=  8; i <<=  8; }
	if (((unsigned int)i) >> 28 == 0) { n +=  4; i <<=  4; }
	if (((unsigned int)i) >> 30 == 0) { n +=  2; i <<=  2; }
	n -= ((unsigned int)i) >> 31;
	return n;
}

inline unsigned char numberOfLeadingZeros_Long(long i) {
	 if (i == 0)
		return 64;
	unsigned char n = 1;
	int x = (int)(((unsigned long)i) >> 32);
	if (x == 0) { n += 32; x = (int)i; }
	if (((unsigned int)x) >> 16 == 0) { n += 16; x <<= 16; }
	if (((unsigned int)x) >> 24 == 0) { n +=  8; x <<=  8; }
	if (((unsigned int)x) >> 28 == 0) { n +=  4; x <<=  4; }
	if (((unsigned int)x) >> 30 == 0) { n +=  2; x <<=  2; }
	n -= ((unsigned int)x) >> 31;
	return n;
}

inline unsigned char getLeadingNumbers_Int(int v1, int v2)
{
	int v = v1 ^ v2;
	return (unsigned char)numberOfLeadingZeros_Int(v);
}

inline unsigned char getLeadingNumbers_Long(long v1, long v2)
{
	long v = v1 ^ v2;
	return (unsigned char)numberOfLeadingZeros_Long(v);
}

/**
 * By default, the endian type is OS endian type.
 * */
inline short bytesToShort(unsigned char* bytes)
{
	lint16 buf;
	memcpy(buf.byte, bytes, 2);
	
	return buf.svalue;
}

inline void shortToBytes(unsigned char* b, short value)
{
	lint16 buf;
	buf.svalue = value;
	memcpy(b, buf.byte, 2);
}

inline int bytesToInt(unsigned char* bytes)
{
	lfloat buf;
	memcpy(buf.byte, bytes, 4);
	return buf.ivalue;
}

inline long bytesToLong(unsigned char* bytes)
{
	ldouble buf;
	memcpy(buf.byte, bytes, 8);
	return buf.lvalue;
}

//the byte to input is in the big-endian format
inline float bytesToFloat(unsigned char* bytes)
{
	lfloat buf;
	memcpy(buf.byte, bytes, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_4bytes(buf.byte);	
	return buf.value;
}

inline void floatToBytes(unsigned char *b, float num)
{
	lfloat buf;
	buf.value = num;
	memcpy(b, buf.byte, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_4bytes(b);		
}

//the byte to input is in the big-endian format
inline double bytesToDouble(unsigned char* bytes)
{
	ldouble buf;
	memcpy(buf.byte, bytes, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_8bytes(buf.byte);
	return buf.value;
}

inline void doubleToBytes(unsigned char *b, double num)
{
	ldouble buf;
	buf.value = num;
	memcpy(b, buf.byte, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		symTransform_8bytes(b);
}


inline int getMaskRightCode(int m) {
	switch (m) {
	case 1:
		return 0x01;
	case 2:
		return 0x03;
	case 3:
		return 0x07;
	case 4:
		return 0x0F;
	case 5:
		return 0x1F;
	case 6:
		return 0x3F;
	case 7:
		return 0X7F;
	case 8:
		return 0XFF;
	default:
		return 0;
	}
}

inline int getLeftMovingCode(int kMod8)
{
	return getMaskRightCode(8 - kMod8);
}

inline int getRightMovingSteps(int kMod8, int resiBitLength) {
	return 8 - kMod8 - resiBitLength;
}

inline int getRightMovingCode(int kMod8, int resiBitLength)
{
	int rightMovingSteps = 8 - kMod8 - resiBitLength;
	if(rightMovingSteps < 0)
	{
		switch(-rightMovingSteps)
		{
		case 1:
			return 0x80;
		case 2:
			return 0xC0;
		case 3:
			return 0xE0;
		case 4:
			return 0xF0;
		case 5:
			return 0xF8;
		case 6:
			return 0xFC;
		case 7:
			return 0XFE;
		default:
			return 0;
		}    		
	}
	else //if(rightMovingSteps >= 0)
	{
		int a = getMaskRightCode(8 - kMod8);
		int b = getMaskRightCode(8 - kMod8 - resiBitLength);
		int c = a - b;
		return c;
	}
}

inline size_t bytesToSize(unsigned char* bytes)
{
	size_t result = bytesToLong_bigEndian(bytes);//8	
	return result;
}

inline void sizeToBytes(unsigned char* outBytes, size_t size)
{
		longToBytes_bigEndian(outBytes, size);//8
}

short* convertByteDataToShortArray(unsigned char* bytes, size_t byteLength);
unsigned short* convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength);
void convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes);
void convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes);
void convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes);
void convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes);
void convertLongArrayToBytes(int64_t* states, size_t stateLength, unsigned char* bytes);
void convertULongArrayToBytes(uint64_t* states, size_t stateLength, unsigned char* bytes);
int computeByteSizePerIntValue(long valueRangeSize);
long computeRangeSize_int(void* oriData, int dataType, size_t size, int64_t* valueRangeSize);
double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue);
float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue);

inline double min_d(double a, double b)
{
	if(a<b)
		return a;
	else
		return b;
}

inline double max_d(double a, double b)
{
	if(a>b)
		return a;
	else
		return b;
}

inline float min_f(float a, float b)
{
	if(a<b)
		return a;
	else
		return b;
}

inline float max_f(float a, float b)
{
	if(a>b)
		return a;
	else
		return b;
}

double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
double getRealPrecision_int(long valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);

inline void compressInt8Value(int8_t tgtValue, int8_t minValue, int byteSize, unsigned char* bytes)
{
	uint8_t data = tgtValue - minValue;
	memcpy(bytes, &data, byteSize); //byteSize==1
}

inline void compressInt16Value(int16_t tgtValue, int16_t minValue, int byteSize, unsigned char* bytes)
{
	uint16_t data = tgtValue - minValue;
	unsigned char tmpBytes[2];
	int16ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 2 - byteSize, byteSize);
}

inline void compressInt32Value(int32_t tgtValue, int32_t minValue, int byteSize, unsigned char* bytes)
{
	uint32_t data = tgtValue - minValue;
	unsigned char tmpBytes[4];
	int32ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 4 - byteSize, byteSize);
}

inline void compressInt64Value(int64_t tgtValue, int64_t minValue, int byteSize, unsigned char* bytes)
{
	uint64_t data = tgtValue - minValue;
	unsigned char tmpBytes[8];
	int64ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 8 - byteSize, byteSize);
}

inline void compressUInt8Value(uint8_t tgtValue, uint8_t minValue, int byteSize, unsigned char* bytes)
{
	uint8_t data = tgtValue - minValue;
	memcpy(bytes, &data, byteSize); //byteSize==1
}

inline void compressUInt16Value(uint16_t tgtValue, uint16_t minValue, int byteSize, unsigned char* bytes)
{
	uint16_t data = tgtValue - minValue;
	unsigned char tmpBytes[2];
	int16ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 2 - byteSize, byteSize);
}

inline void compressUInt32Value(uint32_t tgtValue, uint32_t minValue, int byteSize, unsigned char* bytes)
{
	uint32_t data = tgtValue - minValue;
	unsigned char tmpBytes[4];
	int32ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 4 - byteSize, byteSize);
}

inline void compressUInt64Value(uint64_t tgtValue, uint64_t minValue, int byteSize, unsigned char* bytes)
{
	uint64_t data = tgtValue - minValue;
	unsigned char tmpBytes[8];
	int64ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 8 - byteSize, byteSize);
}

inline int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<8;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}


inline int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<4;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _szp_CompressionToolkit_H  ----- */

