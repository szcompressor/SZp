/**
 *  @file szp_CompressionToolkit.cc
 *  @author Sheng Di
 *  @date Feb, 2025
 *  @brief Compression Toolkit
 *  (C) 2025 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
 
#include <stdlib.h>
#include "szp_defines.h" 	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cstdint>
#include "szp_CompressionToolkit.h"

namespace szp{

short* convertByteDataToShortArray(unsigned char* bytes, size_t byteLength)
{
	szp_lint16 ls;
	size_t i, stateLength = byteLength/2;
	short* states = (short*)malloc(stateLength*sizeof(short));
	if(szp_sysEndianType==szp_dataEndianType)
	{	
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2];
			ls.byte[1] = bytes[i*2+1];
			states[i] = ls.svalue;
		}
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2+1];
			ls.byte[1] = bytes[i*2];
			states[i] = ls.svalue;
		}		
	}
	return states;
} 

unsigned short* convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength)
{
	szp_lint16 ls;
	size_t i, stateLength = byteLength/2;
	unsigned short* states = (unsigned short*)malloc(stateLength*sizeof(unsigned short));
	if(szp_sysEndianType==szp_dataEndianType)
	{	
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2];
			ls.byte[1] = bytes[i*2+1];
			states[i] = ls.usvalue;
		}
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2+1];
			ls.byte[1] = bytes[i*2];
			states[i] = ls.usvalue;
		}		
	}
	return states;
} 

void convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint16 ls;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			ls.svalue = states[i];
			bytes[i*2] = ls.byte[0];
			bytes[i*2+1] = ls.byte[1];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.svalue = states[i];
			bytes[i*2] = ls.byte[1];
			bytes[i*2+1] = ls.byte[0];
		}			
	}
}

void convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint16 ls;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			ls.usvalue = states[i];
			bytes[i*2] = ls.byte[0];
			bytes[i*2+1] = ls.byte[1];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.usvalue = states[i];
			bytes[i*2] = ls.byte[1];
			bytes[i*2+1] = ls.byte[0];
		}			
	}
}

void convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint32 ls;
	size_t index = 0;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.ivalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.ivalue = states[i];
			bytes[index] = ls.byte[3];
			bytes[index+1] = ls.byte[2];
			bytes[index+2] = ls.byte[1];
			bytes[index+3] = ls.byte[0];
		}			
	}
}

void convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint32 ls;
	size_t index = 0;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.uivalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.uivalue = states[i];
			bytes[index] = ls.byte[3];
			bytes[index+1] = ls.byte[2];
			bytes[index+2] = ls.byte[1];
			bytes[index+3] = ls.byte[0];
		}			
	}
}

void convertLongArrayToBytes(int64_t* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint64 ls;
	size_t index = 0;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.lvalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
			bytes[index+4] = ls.byte[4];
			bytes[index+5] = ls.byte[5];
			bytes[index+6] = ls.byte[6];
			bytes[index+7] = ls.byte[7];	
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.lvalue = states[i];
			bytes[index] = ls.byte[7];
			bytes[index+1] = ls.byte[6];
			bytes[index+2] = ls.byte[5];
			bytes[index+3] = ls.byte[4];
			bytes[index+4] = ls.byte[3];
			bytes[index+5] = ls.byte[2];
			bytes[index+6] = ls.byte[1];
			bytes[index+7] = ls.byte[0];	
		}			
	}
}

void convertULongArrayToBytes(uint64_t* states, size_t stateLength, unsigned char* bytes)
{
	szp_lint64 ls;
	size_t index = 0;
	size_t i;
	if(szp_sysEndianType==szp_dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.ulvalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
			bytes[index+4] = ls.byte[4];
			bytes[index+5] = ls.byte[5];
			bytes[index+6] = ls.byte[6];
			bytes[index+7] = ls.byte[7];			
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.ulvalue = states[i];
			bytes[index] = ls.byte[7];
			bytes[index+1] = ls.byte[6];
			bytes[index+2] = ls.byte[5];
			bytes[index+3] = ls.byte[4];
			bytes[index+4] = ls.byte[3];
			bytes[index+5] = ls.byte[2];
			bytes[index+6] = ls.byte[1];
			bytes[index+7] = ls.byte[0];	
		}			
	}
}

int computeByteSizePerIntValue(long valueRangeSize)
{
	if(valueRangeSize<=256)
		return 1;
	else if(valueRangeSize<=65536)
		return 2;
	else if(valueRangeSize<=4294967296) //2^32
		return 4;
	else
		return 8;
}



long computeRangeSize_int(void* oriData, int dataType, size_t size, int64_t* valueRangeSize)
{
	long max = 0, min = 0;

	if(dataType==SZ_UINT8)
	{
		unsigned char* data = static_cast<unsigned char*>(oriData);

		unsigned char min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_INT8)
	{
		char* data = static_cast<char*>(oriData);

		char min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);

	}
	else if(dataType == SZ_UINT16)
	{
		unsigned short* data = static_cast<unsigned short*>(oriData);

		unsigned short min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_INT16)
	{ 
		short* data = static_cast<short*>(oriData);

		short min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_UINT32)
	{
		unsigned int* data = static_cast<unsigned int*>(oriData);

		unsigned int min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_INT32)
	{
		int* data = static_cast<int*>(oriData);

		int min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_UINT64)
	{
		unsigned long* data = static_cast<unsigned long*>(oriData);

		unsigned long min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = static_cast<long>(min_uc);
		max = static_cast<long>(max_uc);
	}
	else if(dataType == SZ_INT64)
	{
		long* data = static_cast<long*>(oriData);

		long min_uc = data[0], max_uc = data[0];
		computeMinMax(data, size, min_uc, max_uc);

		min = min_uc;
		max = max_uc;
	}

	*valueRangeSize = max - min;
	return min;	
}

float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue)
{
	size_t i = 0;
	float min = oriData[0];
	float max = min;
	for(i=1;i<size;i++)
	{
		float data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}

	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue)
{
	size_t i = 0;
	double min = oriData[0];
	double max = min;
	for(i=1;i<size;i++)
	{
		double data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}
	
	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
		precision = relBoundRatio*valueRangeSize;
	else if(errBoundMode==ABS_AND_REL)
		precision = min_d(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==ABS_OR_REL)
		precision = max_d(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==PW_REL)
		precision = 0;
	else
	{
		printf("Error: error-bound-mode is incorrect!\n");
		state = SZ_BERR;
	}
	*status = state;
	return precision;
}

double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
		precision = relBoundRatio*valueRangeSize;
	else if(errBoundMode==ABS_AND_REL)
		precision = min_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==ABS_OR_REL)
		precision = max_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==PW_REL)
		precision = 0;
	else
	{
		printf("Error: error-bound-mode is incorrect!\n");
		state = SZ_BERR;
	}
	*status = state;
	return precision;
}

double getRealPrecision_int(long valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
		precision = relBoundRatio*valueRangeSize;
	else if(errBoundMode==ABS_AND_REL)
		precision = min_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==ABS_OR_REL)
		precision = max_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==PW_REL)
		precision = -1;
	else
	{
		printf("Error: error-bound-mode is incorrect!\n");
		state = SZ_BERR;
	}
	*status = state;
	return precision;
}

}
