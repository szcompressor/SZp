/**
 *  @file szp_TypeManager.c
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>, Sheng Di <sdi1@anl.gov>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include "szp.h"

namespace szp{

size_t Jiajun_save_fixed_length_bits(unsigned int *unsignintArray, size_t intArrayLength, unsigned char *result, unsigned int bit_count)
{

	unsigned int byte_count = 0;
	unsigned int remainder_bit = 0;
	size_t i, j;
	byte_count = bit_count / 8; 
	remainder_bit = bit_count % 8;
	size_t byteLength = 0;
		
	size_t byte_offset = byte_count * intArrayLength;
	if (remainder_bit == 0)
	{
		byteLength = byte_offset;
	}
	else
	{
		byteLength = byte_count * intArrayLength + (remainder_bit * intArrayLength - 1) / 8 + 1;
	}

	size_t n = 0;
	unsigned int tmp;
	i = 0;
	if (byte_count > 0)
	{
		while (n < intArrayLength)
		{
			j = 0;
			tmp = unsignintArray[n];
			tmp >>= remainder_bit;
			while (i < byteLength && j < byte_count)
			{
				result[i] = (unsigned char)tmp;
				i++;
				tmp >>= 8; // right shift by 8 bits to store the next full byte
				j++;
			}
			n++;
		}
	}
	if (remainder_bit > 0)
	{
		if (byte_count > 0)
		{
			for (i = 0; i < intArrayLength; i++)
			{
				unsignintArray[i] <<= (32 - remainder_bit);
				unsignintArray[i] >>= (32 - remainder_bit);
			}
		}
		

		switch (remainder_bit)
		{
		case 1:
			Jiajun_convertUInt2Byte_fast_1b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 2:
			Jiajun_convertUInt2Byte_fast_2b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 3:
			Jiajun_convertUInt2Byte_fast_3b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 4:
			Jiajun_convertUInt2Byte_fast_4b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 5:
			Jiajun_convertUInt2Byte_fast_5b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 6:
			Jiajun_convertUInt2Byte_fast_6b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		case 7:
			Jiajun_convertUInt2Byte_fast_7b_args(unsignintArray, intArrayLength, result + byte_offset);
			break;
		default:
			printf("Error: try to save %d bits\n", remainder_bit);
		}
	}

	return byteLength;
}

size_t convertInt2Byte_fast_1b_args(unsigned char *intArray, size_t intArrayLength, unsigned char *result)
{
	size_t byteLength = 0;
	size_t i, j;
	if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;

	size_t n = 0;
	unsigned int tmp, type;
	for (i = 0; i < byteLength; i++)
	{
		tmp = 0;
		for (j = 0; j < 8 && n < intArrayLength; j++)
		{
			type = intArray[n];
			
			tmp = (tmp | (type << (7 - j)));
			n++;
		}
		result[i] = (unsigned char)tmp;
	}
	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_1b_args(unsigned int *intArray, size_t intArrayLength, unsigned char *result)
{
	size_t byteLength = 0;
	size_t i, j;
	if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;

	size_t n = 0;
	unsigned int tmp, type;
	for (i = 0; i < byteLength; i++)
	{
		tmp = 0;
		for (j = 0; j < 8 && n < intArrayLength; j++)
		{
			type = intArray[n];
			
			tmp = (tmp | (type << (7 - j)));
			n++;
		}
		result[i] = (unsigned char)tmp;
	}
	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_2b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	register unsigned char tmp = 0;
	size_t i, byteLength = 0;
	if (timeStepTypeLength % 4 == 0)
		byteLength = timeStepTypeLength * 2 / 8;
	else
		byteLength = timeStepTypeLength * 2 / 8 + 1;
	size_t n = 0;
	if (timeStepTypeLength % 4 == 0)
	{
		for (i = 0; i < byteLength; i++)
		{
			tmp = 0;

			tmp |= timeStepType[n++] << 6;
			tmp |= timeStepType[n++] << 4;
			tmp |= timeStepType[n++] << 2;
			tmp |= timeStepType[n++];

			

			result[i] = tmp;
		}
	}
	else
	{
		size_t byteLength_ = byteLength - 1;
		for (i = 0; i < byteLength_; i++)
		{
			tmp = 0;

			tmp |= timeStepType[n++] << 6;
			tmp |= timeStepType[n++] << 4;
			tmp |= timeStepType[n++] << 2;
			tmp |= timeStepType[n++];

			

			result[i] = tmp;
		}
		tmp = 0;
		int mod4 = timeStepTypeLength % 4;
		for (int j = 0; j < mod4; j++)
		{
			unsigned char type = timeStepType[n++];
			tmp = tmp | type << (6 - (j << 1));
		}
		result[i] = tmp;
	}

	
	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_3b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 3 / 8;
	else
		byteLength = timeStepTypeLength * 3 / 8 + 1;

	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 5);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 1);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] >> 2);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] << 3);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(result)[i] = tmp;

	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_4b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	size_t i = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 2 == 0)
		byteLength = timeStepTypeLength * 4 / 8;
	else
		byteLength = timeStepTypeLength * 4 / 8 + 1;

	for (n = 0; n < timeStepTypeLength;)
	{
		unsigned char tmp = 0;
		for (int j = 0; j < 2 && n < timeStepTypeLength; j++)
		{
			unsigned int type = timeStepType[n];
			if (j == 0)
				tmp = tmp | (type << 4);
			else 
				tmp = tmp | type;
			n++;
		}
		(result)[i++] = tmp;
	}
	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_5b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 5 / 8;
	else
		byteLength = timeStepTypeLength * 5 / 8 + 1;

	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 3);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 2);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] >> 4);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] >> 1);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] >> 3);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 5);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(result)[i] = tmp;

	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_6b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 6 / 8;
	else
		byteLength = timeStepTypeLength * 6 / 8 + 1;

	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 4;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 4);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 2);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] << 0);
			(result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 3) 
		(result)[i] = tmp;

	return byteLength;
}

size_t Jiajun_convertUInt2Byte_fast_7b_args(unsigned int *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 7 / 8;
	else
		byteLength = timeStepTypeLength * 7 / 8 + 1;

	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 6);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 2);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 5);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 3);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] >> 4);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] >> 3);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 5);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] >> 2);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] >> 1);
			(result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(result)[i] = tmp;

	return byteLength;
}

size_t Jiajun_extract_fixed_length_bits(unsigned char *result, size_t intArrayLength, unsigned int *unsignintArray, unsigned int bit_count)
{
	unsigned int byte_count = 0;
	unsigned int remainder_bit = 0;
	
	size_t i, j;
	byte_count = bit_count / 8; 
	remainder_bit = bit_count % 8;
	size_t byteLength = 0;
	
	
	size_t n = 0;
	unsigned int tmp1, tmp2;
	i = 0;
	size_t byte_offset = byte_count * intArrayLength;

	if (remainder_bit == 0)
	{
		byteLength = byte_offset;
	}
	else
	{
		byteLength = byte_count * intArrayLength + (remainder_bit * intArrayLength - 1) / 8 + 1;
	}

	
	if (remainder_bit > 0)
	{
		switch (remainder_bit)
		{
		case 1:
			Jiajun_convertByte2UInt_fast_1b_args(intArrayLength, result + byte_offset, (intArrayLength - 1) / 8 + 1, unsignintArray);
			break;
		case 2:
			Jiajun_convertByte2UInt_fast_2b_args(intArrayLength, result + byte_offset, (intArrayLength * 2 - 1) / 8 + 1, unsignintArray);
			break;
		case 3:
			Jiajun_convertByte2UInt_fast_3b_args(intArrayLength, result + byte_offset, (intArrayLength * 3 - 1) / 8 + 1, unsignintArray);
			break;
		case 4:
			Jiajun_convertByte2UInt_fast_4b_args(intArrayLength, result + byte_offset, (intArrayLength * 4 - 1) / 8 + 1, unsignintArray);
			break;
		case 5:
			Jiajun_convertByte2UInt_fast_5b_args(intArrayLength, result + byte_offset, (intArrayLength * 5 - 1) / 8 + 1, unsignintArray);
			break;
		case 6:
			Jiajun_convertByte2UInt_fast_6b_args(intArrayLength, result + byte_offset, (intArrayLength * 6 - 1) / 8 + 1, unsignintArray);
			break;
		case 7:
			Jiajun_convertByte2UInt_fast_7b_args(intArrayLength, result + byte_offset, (intArrayLength * 7 - 1) / 8 + 1, unsignintArray);
			break;
		default:
			printf("Error: try to extract %d bits\n", remainder_bit);
		}
	}
	if (byte_count > 0)
	{
		if(remainder_bit == 0)
		{
			memset(unsignintArray, 0 , intArrayLength * sizeof(unsigned int));
		}
		while (i < intArrayLength)
		{
			j = 0;
			tmp1 = 0;
			tmp2 = 0;
			while (j < byte_count && n < byteLength)
			{
				tmp1 = result[n];
				tmp1 <<= (8 * j);
				tmp2 = (tmp2 | tmp1);
				n++;
				j++;
			}
			tmp2 <<= remainder_bit; 
			unsignintArray[i] = (unsignintArray[i] | tmp2);
			i++;
		}
	}

	return byteLength;
}

void Jiajun_convertByte2UInt_fast_1b_args(size_t intArrayLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t n = 0, i;
	unsigned int tmp;
	for (i = 0; i < byteArrayLength - 1; i++)
	{
		tmp = byteArray[i];
		intArray[n++] = (tmp & 0x80) >> 7;
		intArray[n++] = (tmp & 0x40) >> 6;
		intArray[n++] = (tmp & 0x20) >> 5;
		intArray[n++] = (tmp & 0x10) >> 4;
		intArray[n++] = (tmp & 0x08) >> 3;
		intArray[n++] = (tmp & 0x04) >> 2;
		intArray[n++] = (tmp & 0x02) >> 1;
		intArray[n++] = (tmp & 0x01) >> 0;
	}

	tmp = byteArray[i];
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x80) >> 7;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x40) >> 6;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x20) >> 5;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x10) >> 4;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x08) >> 3;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x04) >> 2;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x02) >> 1;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x01) >> 0;
}

void Jiajun_convertByte2UInt_fast_2b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{

	size_t i, n = 0;

	int mod4 = stepLength % 4;
	if (mod4 == 0)
	{
		for (i = 0; i < byteArrayLength; i++)
		{
			unsigned char tmp = byteArray[i];
			intArray[n++] = (tmp & 0xC0) >> 6;
			intArray[n++] = (tmp & 0x30) >> 4;
			intArray[n++] = (tmp & 0x0C) >> 2;
			intArray[n++] = tmp & 0x03;
		}
	}
	else
	{
		size_t t = byteArrayLength - 1;
		for (i = 0; i < t; i++)
		{
			unsigned char tmp = byteArray[i];
			intArray[n++] = (tmp & 0xC0) >> 6;
			intArray[n++] = (tmp & 0x30) >> 4;
			intArray[n++] = (tmp & 0x0C) >> 2;
			intArray[n++] = tmp & 0x03;
		}
		unsigned char tmp = byteArray[i];
		switch (mod4)
		{
		case 1:
			intArray[n++] = (tmp & 0xC0) >> 6;
			break;
		case 2:
			intArray[n++] = (tmp & 0xC0) >> 6;
			intArray[n++] = (tmp & 0x30) >> 4;
			break;
		case 3:
			intArray[n++] = (tmp & 0xC0) >> 6;
			intArray[n++] = (tmp & 0x30) >> 4;
			intArray[n++] = (tmp & 0x0C) >> 2;
			break;
		}
	}
}

void Jiajun_convertByte2UInt_fast_3b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			intArray[n++] = (tmp & 0xE0) >> 5;
			break;
		case 1:
			intArray[n++] = (tmp & 0x1C) >> 2;
			break;
		case 2:
			ii = (tmp & 0x03) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			intArray[n++] = ii;
			break;
		case 3:
			intArray[n++] = (tmp & 0x70) >> 4;
			break;
		case 4:
			intArray[n++] = (tmp & 0x0E) >> 1;
			break;
		case 5:
			ii = (tmp & 0x01) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			intArray[n++] = ii;
			break;
		case 6:
			intArray[n++] = (tmp & 0x38) >> 3;
			break;
		case 7:
			intArray[n++] = (tmp & 0x07);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

void Jiajun_convertByte2UInt_fast_4b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t i = 0, n = 0;
	unsigned char tmp;
	for (i = 0; i < byteArrayLength; i++)
	{
		tmp = byteArray[i];
		intArray[n++] = (tmp & 0xF0) >> 4;
		if (n == stepLength)
			break;
		intArray[n++] = tmp & 0x0F;
	}
}

void Jiajun_convertByte2UInt_fast_5b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			intArray[n++] = (tmp & 0xF8) >> 3;
			break;
		case 1:
			ii = (tmp & 0x07) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			intArray[n++] = ii;
			break;
		case 2:
			intArray[n++] = (tmp & 0x3E) >> 1;
			break;
		case 3:
			ii = (tmp & 0x01) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			intArray[n++] = ii;
			break;
		case 4:
			ii = (tmp & 0x0F) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			intArray[n++] = ii;
			break;
		case 5:
			intArray[n++] = (tmp & 0x7C) >> 2;
			break;
		case 6:
			ii = (tmp & 0x03) << 3;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xE0) >> 5;
			intArray[n++] = ii;
			break;
		case 7:
			intArray[n++] = (tmp & 0x1F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

void Jiajun_convertByte2UInt_fast_6b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 4)
		{
		case 0:
			intArray[n++] = (tmp & 0xFC) >> 2;
			break;
		case 1:
			ii = (tmp & 0x03) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			intArray[n++] = ii;
			break;
		case 2:
			ii = (tmp & 0x0F) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			intArray[n++] = ii;
			break;
		case 3:
			intArray[n++] = (tmp & 0x3F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

void Jiajun_convertByte2UInt_fast_7b_args(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned int *intArray)
{
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			intArray[n++] = (tmp & 0xFE) >> 1;
			break;
		case 1:
			ii = (tmp & 0x01) << 6;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xFC) >> 2;
			intArray[n++] = ii;
			break;
		case 2:
			ii = (tmp & 0x03) << 5;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF8) >> 3;
			intArray[n++] = ii;
			break;
		case 3:
			ii = (tmp & 0x07) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			intArray[n++] = ii;
			break;
		case 4:
			ii = (tmp & 0x0F) << 3;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xE0) >> 5;
			intArray[n++] = ii;
			break;
		case 5:
			ii = (tmp & 0x1F) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			intArray[n++] = ii;
			break;
		case 6:
			ii = (tmp & 0x3F) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			intArray[n++] = ii;
			break;
		case 7:
			intArray[n++] = (tmp & 0x7F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

size_t convertIntArray2ByteArray_fast_1b_args(unsigned char *intArray, size_t intArrayLength, unsigned char *result)
{
	size_t byteLength = 0;
	size_t i, j;
	if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;

	size_t n = 0;
	int tmp, type;
	for (i = 0; i < byteLength; i++)
	{
		tmp = 0;
		for (j = 0; j < 8 && n < intArrayLength; j++)
		{
			type = intArray[n];
			
			tmp = (tmp | (type << (7 - j)));
			n++;
		}
		result[i] = (unsigned char)tmp;
	}
	return byteLength;
}

size_t convertIntArray2ByteArray_fast_1b(unsigned char *intArray, size_t intArrayLength, unsigned char **result)
{
	size_t byteLength = 0;
	size_t i, j;
	if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	size_t n = 0;
	int tmp, type;
	for (i = 0; i < byteLength; i++)
	{
		tmp = 0;
		for (j = 0; j < 8 && n < intArrayLength; j++)
		{
			type = intArray[n];
			if (type == 1)
				tmp = (tmp | (1 << (7 - j)));
			n++;
		}
		(*result)[i] = (unsigned char)tmp;
	}
	return byteLength;
}

size_t convertIntArray2ByteArray_fast_1b_to_result(unsigned char *intArray, size_t intArrayLength, unsigned char *result)
{
	size_t byteLength = 0;
	size_t i, j;
	if (intArrayLength % 8 == 0)
		byteLength = intArrayLength / 8;
	else
		byteLength = intArrayLength / 8 + 1;

	size_t n = 0;
	int tmp, type;
	for (i = 0; i < byteLength; i++)
	{
		tmp = 0;
		for (j = 0; j < 8 && n < intArrayLength; j++)
		{
			type = intArray[n];
			if (type == 1)
				tmp = (tmp | (1 << (7 - j)));
			n++;
		}
		result[i] = (unsigned char)tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_1b_args(size_t intArrayLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char *intArray)
{
	size_t n = 0, i;
	unsigned int tmp;
	for (i = 0; i < byteArrayLength - 1; i++)
	{
		tmp = byteArray[i];
		intArray[n++] = (tmp & 0x80) >> 7;
		intArray[n++] = (tmp & 0x40) >> 6;
		intArray[n++] = (tmp & 0x20) >> 5;
		intArray[n++] = (tmp & 0x10) >> 4;
		intArray[n++] = (tmp & 0x08) >> 3;
		intArray[n++] = (tmp & 0x04) >> 2;
		intArray[n++] = (tmp & 0x02) >> 1;
		intArray[n++] = (tmp & 0x01) >> 0;
	}

	tmp = byteArray[i];
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x80) >> 7;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x40) >> 6;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x20) >> 5;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x10) >> 4;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x08) >> 3;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x04) >> 2;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x02) >> 1;
	if (n == intArrayLength)
		return;
	intArray[n++] = (tmp & 0x01) >> 0;
}

void convertByteArray2IntArray_fast_1b(size_t intArrayLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (intArrayLength > byteArrayLength * 8)
	{
		printf("Error: intArrayLength > byteArrayLength*8\n");
		printf("intArrayLength=%zu, byteArrayLength = %zu", intArrayLength, byteArrayLength);
		exit(0);
	}
	if (intArrayLength > 0)
		*intArray = (unsigned char *)malloc(intArrayLength * sizeof(unsigned char));
	else
		*intArray = NULL;

	size_t n = 0, i;
	int tmp;
	for (i = 0; i < byteArrayLength - 1; i++)
	{
		tmp = byteArray[i];
		(*intArray)[n++] = (tmp & 0x80) >> 7;
		(*intArray)[n++] = (tmp & 0x40) >> 6;
		(*intArray)[n++] = (tmp & 0x20) >> 5;
		(*intArray)[n++] = (tmp & 0x10) >> 4;
		(*intArray)[n++] = (tmp & 0x08) >> 3;
		(*intArray)[n++] = (tmp & 0x04) >> 2;
		(*intArray)[n++] = (tmp & 0x02) >> 1;
		(*intArray)[n++] = (tmp & 0x01) >> 0;
	}

	tmp = byteArray[i];
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x80) >> 7;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x40) >> 6;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x20) >> 5;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x10) >> 4;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x08) >> 3;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x04) >> 2;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x02) >> 1;
	if (n == intArrayLength)
		return;
	(*intArray)[n++] = (tmp & 0x01) >> 0;
}

inline size_t convertIntArray2ByteArray_fast_2b_args(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char *result)
{
	register unsigned char tmp = 0;
	size_t i, byteLength = 0;
	if (timeStepTypeLength % 4 == 0)
		byteLength = timeStepTypeLength * 2 / 8;
	else
		byteLength = timeStepTypeLength * 2 / 8 + 1;
	size_t n = 0;
	if (timeStepTypeLength % 4 == 0)
	{
		for (i = 0; i < byteLength; i++)
		{
			tmp = 0;

			tmp |= timeStepType[n++] << 6;
			tmp |= timeStepType[n++] << 4;
			tmp |= timeStepType[n++] << 2;
			tmp |= timeStepType[n++];

			

			result[i] = tmp;
		}
	}
	else
	{
		size_t byteLength_ = byteLength - 1;
		for (i = 0; i < byteLength_; i++)
		{
			tmp = 0;

			tmp |= timeStepType[n++] << 6;
			tmp |= timeStepType[n++] << 4;
			tmp |= timeStepType[n++] << 2;
			tmp |= timeStepType[n++];

			

			result[i] = tmp;
		}
		tmp = 0;
		int mod4 = timeStepTypeLength % 4;
		for (int j = 0; j < mod4; j++)
		{
			unsigned char type = timeStepType[n++];
			tmp = tmp | type << (6 - (j << 1));
		}
		result[i] = tmp;
	}

	
	return byteLength;
}

/**
 * little endian
 * [01|10|11|00|....]-->[01|10|11|00][....]
 * @param timeStepType
 * @return
 */
size_t convertIntArray2ByteArray_fast_2b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i, j, byteLength = 0;
	if (timeStepTypeLength % 4 == 0)
		byteLength = timeStepTypeLength * 2 / 8;
	else
		byteLength = timeStepTypeLength * 2 / 8 + 1;
	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	size_t n = 0;
	for (i = 0; i < byteLength; i++)
	{
		int tmp = 0;
		for (j = 0; j < 4 && n < timeStepTypeLength; j++)
		{
			int type = timeStepType[n];
			switch (type)
			{
			case 0:

				break;
			case 1:
				tmp = (tmp | (1 << (6 - j * 2)));
				break;
			case 2:
				tmp = (tmp | (2 << (6 - j * 2)));
				break;
			case 3:
				tmp = (tmp | (3 << (6 - j * 2)));
				break;
			default:
				printf("Error: wrong timestep type...: type[%zu]=%d\n", n, type);
				exit(0);
			}
			n++;
		}
		(*result)[i] = (unsigned char)tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 4)
	{
		printf("Error: stepLength > byteArray.length*4\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i, n = 0;

	int mod4 = stepLength % 4;
	if (mod4 == 0)
	{
		for (i = 0; i < byteArrayLength; i++)
		{
			unsigned char tmp = byteArray[i];
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;
			(*intArray)[n++] = tmp & 0x03;
		}
	}
	else
	{
		size_t t = byteArrayLength - 1;
		for (i = 0; i < t; i++)
		{
			unsigned char tmp = byteArray[i];
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;
			(*intArray)[n++] = tmp & 0x03;
		}
		unsigned char tmp = byteArray[i];
		switch (mod4)
		{
		case 1:
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			break;
		case 2:
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			break;
		case 3:
			(*intArray)[n++] = (tmp & 0xC0) >> 6;
			(*intArray)[n++] = (tmp & 0x30) >> 4;
			(*intArray)[n++] = (tmp & 0x0C) >> 2;
			break;
		}
	}
}

size_t convertIntArray2ByteArray_fast_3b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 3 / 8;
	else
		byteLength = timeStepTypeLength * 3 / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 5);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 1);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] >> 2);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] << 3);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(*result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(*result)[i] = tmp;

	return byteLength;
}

void convertByteArray2IntArray_fast_3b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 8 / 3)
	{
		printf("Error: stepLength > byteArray.length*8/3, impossible case unless bugs elsewhere.\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			(*intArray)[n++] = (tmp & 0xE0) >> 5;
			break;
		case 1:
			(*intArray)[n++] = (tmp & 0x1C) >> 2;
			break;
		case 2:
			ii = (tmp & 0x03) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			(*intArray)[n++] = ii;
			break;
		case 3:
			(*intArray)[n++] = (tmp & 0x70) >> 4;
			break;
		case 4:
			(*intArray)[n++] = (tmp & 0x0E) >> 1;
			break;
		case 5:
			ii = (tmp & 0x01) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			(*intArray)[n++] = ii;
			break;
		case 6:
			(*intArray)[n++] = (tmp & 0x38) >> 3;
			break;
		case 7:
			(*intArray)[n++] = (tmp & 0x07);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

size_t convertIntArray2ByteArray_fast_4b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 2 == 0)
		byteLength = timeStepTypeLength * 4 / 8;
	else
		byteLength = timeStepTypeLength * 4 / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		unsigned char tmp = 0;
		for (int j = 0; j < 2 && n < timeStepTypeLength; j++)
		{
			int type = timeStepType[n];
			if (j == 0)
				tmp = tmp | (type << 4);
			else 
				tmp = tmp | type;
			n++;
		}
		(*result)[i++] = tmp;
	}
	return byteLength;
}

void convertByteArray2IntArray_fast_4b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 2)
	{
		printf("Error: stepLength > byteArray.length*2, impossible case unless bugs elsewhere.\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i = 0, n = 0;
	unsigned char tmp;
	for (i = 0; i < byteArrayLength; i++)
	{
		tmp = byteArray[i];
		(*intArray)[n++] = (tmp & 0xF0) >> 4;
		if (n == stepLength)
			break;
		(*intArray)[n++] = tmp & 0x0F;
	}
}

size_t convertIntArray2ByteArray_fast_5b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 5 / 8;
	else
		byteLength = timeStepTypeLength * 5 / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 3);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 2);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] >> 4);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] >> 1);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] >> 3);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 5);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(*result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(*result)[i] = tmp;

	return byteLength;
}

void convertByteArray2IntArray_fast_5b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 8 / 5)
	{
		printf("Error: stepLength > byteArray.length*8/5, impossible case unless bugs elsewhere.\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			(*intArray)[n++] = (tmp & 0xF8) >> 3;
			break;
		case 1:
			ii = (tmp & 0x07) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			(*intArray)[n++] = ii;
			break;
		case 2:
			(*intArray)[n++] = (tmp & 0x3E) >> 1;
			break;
		case 3:
			ii = (tmp & 0x01) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			(*intArray)[n++] = ii;
			break;
		case 4:
			ii = (tmp & 0x0F) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			(*intArray)[n++] = ii;
			break;
		case 5:
			(*intArray)[n++] = (tmp & 0x7C) >> 2;
			break;
		case 6:
			ii = (tmp & 0x03) << 3;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xE0) >> 5;
			(*intArray)[n++] = ii;
			break;
		case 7:
			(*intArray)[n++] = (tmp & 0x1F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

size_t convertIntArray2ByteArray_fast_6b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 6 / 8;
	else
		byteLength = timeStepTypeLength * 6 / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 4;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 2);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 4);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 2);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] << 0);
			(*result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 3) 
		(*result)[i] = tmp;

	return byteLength;
}

void convertByteArray2IntArray_fast_6b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 8 / 6)
	{
		printf("Error: stepLength > byteArray.length*8/6, impossible case unless bugs elsewhere.\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 4)
		{
		case 0:
			(*intArray)[n++] = (tmp & 0xFC) >> 2;
			break;
		case 1:
			ii = (tmp & 0x03) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			(*intArray)[n++] = ii;
			break;
		case 2:
			ii = (tmp & 0x0F) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			(*intArray)[n++] = ii;
			break;
		case 3:
			(*intArray)[n++] = (tmp & 0x3F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

size_t convertIntArray2ByteArray_fast_7b(unsigned char *timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
	size_t i = 0, k = 0, byteLength = 0, n = 0;
	if (timeStepTypeLength % 8 == 0)
		byteLength = timeStepTypeLength * 7 / 8;
	else
		byteLength = timeStepTypeLength * 7 / 8 + 1;

	if (byteLength > 0)
		*result = (unsigned char *)malloc(byteLength * sizeof(unsigned char));
	else
		*result = NULL;
	unsigned char tmp = 0;
	for (n = 0; n < timeStepTypeLength; n++)
	{
		k = n % 8;
		switch (k)
		{
		case 0:
			tmp = tmp | (timeStepType[n] << 1);
			break;
		case 1:
			tmp = tmp | (timeStepType[n] >> 6);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 2);
			break;
		case 2:
			tmp = tmp | (timeStepType[n] >> 5);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 3);
			break;
		case 3:
			tmp = tmp | (timeStepType[n] >> 4);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 4);
			break;
		case 4:
			tmp = tmp | (timeStepType[n] >> 3);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 5);
			break;
		case 5:
			tmp = tmp | (timeStepType[n] >> 2);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 6);
			break;
		case 6:
			tmp = tmp | (timeStepType[n] >> 1);
			(*result)[i++] = tmp;
			tmp = 0 | (timeStepType[n] << 7);
			break;
		case 7:
			tmp = tmp | (timeStepType[n] << 0);
			(*result)[i++] = tmp;
			tmp = 0;
			break;
		}
	}
	if (k != 7) 
		(*result)[i] = tmp;

	return byteLength;
}

void convertByteArray2IntArray_fast_7b(size_t stepLength, unsigned char *byteArray, size_t byteArrayLength, unsigned char **intArray)
{
	if (stepLength > byteArrayLength * 8 / 7)
	{
		printf("Error: stepLength > byteArray.length*8/7, impossible case unless bugs elsewhere.\n");
		printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
		exit(0);
	}
	if (stepLength > 0)
		*intArray = (unsigned char *)malloc(stepLength * sizeof(unsigned char));
	else
		*intArray = NULL;
	size_t i = 0, ii = 0, n = 0;
	unsigned char tmp = byteArray[i];
	for (n = 0; n < stepLength;)
	{
		switch (n % 8)
		{
		case 0:
			(*intArray)[n++] = (tmp & 0xFE) >> 1;
			break;
		case 1:
			ii = (tmp & 0x01) << 6;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xFC) >> 2;
			(*intArray)[n++] = ii;
			break;
		case 2:
			ii = (tmp & 0x03) << 5;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF8) >> 3;
			(*intArray)[n++] = ii;
			break;
		case 3:
			ii = (tmp & 0x07) << 4;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xF0) >> 4;
			(*intArray)[n++] = ii;
			break;
		case 4:
			ii = (tmp & 0x0F) << 3;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xE0) >> 5;
			(*intArray)[n++] = ii;
			break;
		case 5:
			ii = (tmp & 0x1F) << 2;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0xC0) >> 6;
			(*intArray)[n++] = ii;
			break;
		case 6:
			ii = (tmp & 0x3F) << 1;
			i++;
			tmp = byteArray[i];
			ii |= (tmp & 0x80) >> 7;
			(*intArray)[n++] = ii;
			break;
		case 7:
			(*intArray)[n++] = (tmp & 0x7F);
			i++;
			tmp = byteArray[i];
			break;
		}
	}
}

inline int getLeftMovingSteps(size_t k, unsigned char resiBitLength)
{
	return 8 - k % 8 - resiBitLength;
}

}
