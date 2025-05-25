/**
 *  @file szp_rw.c
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief io interface for fortrance
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>

#include "szp_rw.h"
#include "szp.h"
#include "szp_CompressionToolkit.h"

using namespace szp;

int szp_checkFileExistance(char* filePath)
{
	if( access( filePath, F_OK ) != -1 ) {
		// file exists
		return 1;
	} else {
		// file doesn't exist
		return 0;
	}	
}

float** szp_create2DArray_float(size_t m, size_t n)
{
	size_t i=0;
	float **data = (float**)malloc(sizeof(float*)*m);
	for(i=0;i<m;i++)
		data[i] = (float*)malloc(sizeof(float)*n);
	return data;
}

void szp_free2DArray_float(float** data, size_t m)
{
	size_t i = 0;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

float*** szp_create3DArray_float(size_t p, size_t m, size_t n)
{
	size_t i = 0, j = 0;
	float ***data = (float***)malloc(sizeof(float**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (float**)malloc(sizeof(float*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (float*)malloc(sizeof(float)*n);
	}
	return data;
}

void szp_free3DArray_float(float*** data, size_t p, size_t m)
{
	size_t i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

double** szp_create2DArray_double(size_t m, size_t n)
{
	size_t i=0;
	double **data = (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
			data[i] = (double*)malloc(sizeof(double)*n);
			
	return data;
}

void szp_free2DArray_double(double** data, size_t m)
{
	size_t i;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

double*** szp_create3DArray_double(size_t p, size_t m, size_t n)
{
	size_t i = 0, j = 0;
	double ***data = (double***)malloc(sizeof(double**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (double**)malloc(sizeof(double*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (double*)malloc(sizeof(double)*n);
	}
	return data;
}

void szp_free3DArray_double(double*** data, size_t p, size_t m)
{
	size_t i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

size_t szp_checkFileSize(char *srcFilePath, int *status)
{
	size_t filesize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return -1;
	}
	fseek(pFile, 0, SEEK_END);
    filesize = ftell(pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return filesize;
}

unsigned char *szp_readByteData(char *srcFilePath, size_t *byteLength, int *status)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return 0;
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = ftell(pFile);
    fclose(pFile);
    
    unsigned char *byteBuf = ( unsigned char *)malloc((*byteLength)*sizeof(unsigned char)); //sizeof(char)==1
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return 0;
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return byteBuf;
}

double *szp_readDoubleData(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		double *daBuf = szp_readDoubleData_systemEndian(srcFilePath, nbEle,&state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state==SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		double *daBuf = (double *)malloc(byteLength);
		*nbEle = byteLength/8;
		
		ldouble buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*8;
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}


int8_t *szp_readInt8Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	int8_t *daBuf = szp_readInt8Data_systemEndian(srcFilePath, nbEle, &state);
	*status = state;
	return daBuf;
}

int16_t *szp_readInt16Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		int16_t *daBuf = szp_readInt16Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		int16_t *daBuf = (int16_t *)malloc(byteLength);
		*nbEle = byteLength/2;

		lint16 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 1;//*2
			memcpy(buf.byte, bytes+j, 2);
			symTransform_2bytes(buf.byte);
			daBuf[i] = buf.svalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint16_t *szp_readUInt16Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint16_t *daBuf = szp_readUInt16Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		uint16_t *daBuf = (uint16_t *)malloc(byteLength);
		*nbEle = byteLength/2;

		lint16 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 1;//*2
			memcpy(buf.byte, bytes+j, 2);
			symTransform_2bytes(buf.byte);
			daBuf[i] = buf.usvalue;
		}
		free(bytes);
		return daBuf;
	}
}

int32_t *szp_readInt32Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		int32_t *daBuf = szp_readInt32Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		int32_t *daBuf = (int32_t *)malloc(byteLength);
		*nbEle = byteLength/4;

		lint32 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.ivalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint32_t *szp_readUInt32Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint32_t *daBuf = szp_readUInt32Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		uint32_t *daBuf = (uint32_t *)malloc(byteLength);
		*nbEle = byteLength/4;

		lint32 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 2; //*4
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.uivalue;
		}
		free(bytes);
		return daBuf;
	}
}

int64_t *szp_readInt64Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		int64_t *daBuf = szp_readInt64Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		int64_t *daBuf = (int64_t *)malloc(byteLength);
		*nbEle = byteLength/8;

		lint64 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 3; //*8
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.lvalue;
		}
		free(bytes);
		return daBuf;
	}
}

uint64_t *szp_readUInt64Data(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		uint64_t *daBuf = szp_readUInt64Data_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		uint64_t *daBuf = (uint64_t *)malloc(byteLength);
		*nbEle = byteLength/8;

		lint64 buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i << 3; //*8
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.ulvalue;
		}
		free(bytes);
		return daBuf;
	}
}


float *szp_readFloatData(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = szp_readFloatData_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;
		
		size_t byteLength;
		unsigned char* bytes = szp_readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;
		
		lfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

double *szp_readDoubleData_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = inSize/8; //only support double in this version
    fclose(pFile);
    
    double *daBuf = (double *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return NULL;
    }
    fread(daBuf, 8, *nbEle, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return daBuf;
}


int8_t *szp_readInt8Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize;
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	int8_t *daBuf = (int8_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 1, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;
}


int16_t *szp_readInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/2; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	int16_t *daBuf = (int16_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 2, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;	
}

uint16_t *szp_readUInt16Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/2; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	uint16_t *daBuf = (uint16_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 2, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;	
}

int32_t *szp_readInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/4; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	int32_t *daBuf = (int32_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 4, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;	
}

uint32_t *szp_readUInt32Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/4; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	uint32_t *daBuf = (uint32_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 4, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;	
}

int64_t *szp_readInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/8; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	int64_t *daBuf = (int64_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 8, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;
}

uint64_t *szp_readUInt64Data_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return NULL;
	}
	fseek(pFile, 0, SEEK_END);
	inSize = ftell(pFile);
	*nbEle = inSize/8; 
	fclose(pFile);

	if(inSize<=0)
	{
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}

	uint64_t *daBuf = (uint64_t *)malloc(inSize);

	pFile = fopen(srcFilePath, "rb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 2\n");
		*status = SZ_FERR;
		return NULL;
	}
	fread(daBuf, 8, *nbEle, pFile);
	fclose(pFile);
	*status = SZ_SCES;
	return daBuf;
}

float *szp_readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = inSize/4; 
    fclose(pFile);
    
    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}
    
    float *daBuf = (float *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return NULL;
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return daBuf;
}

void szp_writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
    
    //printf("DEBUG: About to write %zu bytes to %s (bytes=%p)\n", byteLength, tgtFilePath, bytes);
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = SZ_SCES;
}

void szp_writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status)
{
	size_t i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
    
    for(i = 0;i<nbEle;i++)
	{
		sprintf(s,"%.20G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = SZ_SCES;
}

void szp_writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status)
{
	size_t i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
   
    for(i = 0;i<nbEle;i++)
	{
		//printf("i=%d\n",i);
		//printf("data[i]=%f\n",data[i]);
		sprintf(s,"%.30G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = SZ_SCES;
}

void szp_writeData(void *data, int dataType, size_t nbEle, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	if(dataType == SZ_FLOAT)
	{
		float* dataArray = (float *)data;
		szp_writeFloatData(dataArray, nbEle, tgtFilePath, &state);
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* dataArray = (double *)data;
		szp_writeDoubleData(dataArray, nbEle, tgtFilePath, &state);	
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		*status = SZ_TERR; //wrong type
		return;
	}
	*status = state;
}

void szp_writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status)
{
	size_t i = 0; 
	int state = SZ_SCES;
	lfloat buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
	for(i=0;i<nbEle;i++)
	{
		buf.value = data[i];
		bytes[i*4+0] = buf.byte[0];
		bytes[i*4+1] = buf.byte[1];
		bytes[i*4+2] = buf.byte[2];
		bytes[i*4+3] = buf.byte[3];					
	}

	size_t byteLength = nbEle*sizeof(float);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status)
{
	size_t i = 0, index = 0; 
	int state = SZ_SCES;
	ldouble buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
	for(i=0;i<nbEle;i++)
	{
		index = i*8;
		buf.value = data[i];
		bytes[index+0] = buf.byte[0];
		bytes[index+1] = buf.byte[1];
		bytes[index+2] = buf.byte[2];
		bytes[index+3] = buf.byte[3];
		bytes[index+4] = buf.byte[4];
		bytes[index+5] = buf.byte[5];
		bytes[index+6] = buf.byte[6];
		bytes[index+7] = buf.byte[7];
	}

	size_t byteLength = nbEle*sizeof(double);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeShortData_inBytes(short *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertShortArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeUShortData_inBytes(unsigned short *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertUShortArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeIntData_inBytes(int *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*4;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertIntArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeUIntData_inBytes(unsigned int *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*4;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertUIntArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeLongData_inBytes(int64_t *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*8;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertLongArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void szp_writeULongData_inBytes(uint64_t *states, size_t stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	size_t byteLength = stateLength*8;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertULongArrayToBytes(states, stateLength, bytes);
	szp_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

unsigned short* szp_readShortData(char *srcFilePath, size_t *dataLength, int *status)
{
	size_t byteLength = 0; 
	int state = SZ_SCES;
	unsigned char * bytes = szp_readByteData(srcFilePath, &byteLength, &state);
	*dataLength = byteLength/2;
	unsigned short* states = convertByteDataToUShortArray(bytes, byteLength);
	free(bytes);
	*status = state;
	return states;
}

void szp_writeStrings(int nbStr, char *str[], char *tgtFilePath, int *status)
{
	int i = 0;
	char s[256];
	FILE *pFile = fopen(tgtFilePath, "wb");
	if (pFile == NULL)
	{
		printf("Failed to open input file. 3\n");
		*status = SZ_FERR;
		return;
	}

	for(i = 0;i<nbStr;i++)
	{
		sprintf(s,"%s\n",str[i]);
		fputs(s, pFile);
	}

	fclose(pFile);
	*status = SZ_SCES;
}
