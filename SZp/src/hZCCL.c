/**
 *  @file hZCCL.c
 *  @author Sheng Di
 *  @date Jan, 2022
 *  @brief 
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "hZCCL.h"
#include "hZCCL_rw.h"

int versionNumber[4] = {hZCCL_VER_MAJOR,hZCCL_VER_MINOR,hZCCL_VER_BUILD,hZCCL_VER_REVISION};

int dataEndianType = LITTLE_ENDIAN_DATA; //*endian type of the data read from disk
int sysEndianType = LITTLE_ENDIAN_SYSTEM; //*sysEndianType is actually set automatically.

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int dimension;
	if(r1==0)
	{
		dimension = 0;
	}
	else if(r2==0)
	{
		dimension = 1;
	}
	else if(r3==0)
	{
		dimension = 2;
	}
	else if(r4==0)
	{
		dimension = 3;
	}
	else if(r5==0)
	{
		dimension = 4;
	}
	else
	{
		dimension = 5;
	}
	return dimension;
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0)
	{
		dataLength = 0;
	}
	else if(r2==0)
	{
		dataLength = r1;
	}
	else if(r3==0)
	{
		dataLength = r1*r2;
	}
	else if(r4==0)
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0)
	{
		dataLength = r1*r2*r3*r4;
	}
	else
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/**
 * @brief		check dimension and correct it if needed
 * @return 	0 (didn't change dimension)
 * 					1 (dimension is changed)
 * 					2 (dimension is problematic)
 **/
int filterDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* correctedDimension)
{
	int dimensionCorrected = 0;
	int dim = computeDimension(r5, r4, r3, r2, r1);
	correctedDimension[0] = r1;
	correctedDimension[1] = r2;
	correctedDimension[2] = r3;
	correctedDimension[3] = r4;
	correctedDimension[4] = r5;
	size_t* c = correctedDimension;
	if(dim==1)
	{
		if(r1<1)
			return 2;
	}
	else if(dim==2)
	{
		if(r2==1)
		{
			c[1]= 0;
			dimensionCorrected = 1;
		}	
		if(r1==1) //remove this dimension
		{
			c[0] = c[1]; 
			c[1] = c[2];
			dimensionCorrected = 1;
		}
	}
	else if(dim==3)
	{
		if(r3==1)
		{
			c[2] = 0;
			dimensionCorrected = 1;
		}	
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			dimensionCorrected = 1;
		}
	}
	else if(dim==4)
	{
		if(r4==1)
		{
			c[3] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			dimensionCorrected = 1;
		}
	}
	else if(dim==5)
	{
		if(r5==1)
		{
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r4==1)
		{
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r3==1)
		{
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r2==1)
		{
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
		if(r1==1)
		{
			c[0] = c[1];
			c[1] = c[2];
			c[2] = c[3];
			c[3] = c[4];
			c[4] = 0;
			dimensionCorrected = 1;
		}
	}
	
	return dimensionCorrected;
	
}

