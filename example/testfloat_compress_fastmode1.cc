/**
 *  @file szp_Float.h
 *  @author Jiajun Huang, Sheng Di
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include "szp.h"
#ifdef _OPENMP
#include "omp.h"
#endif
struct timeval startTime;
struct timeval endTime;   /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

void cost_start()
{
  totalCost = 0;
  gettimeofday(&costStart, NULL);
}

void cost_end()
{
  double elapsed;
  struct timeval costEnd;
  gettimeofday(&costEnd, NULL);
  elapsed = ((costEnd.tv_sec * 1000000 + costEnd.tv_usec) - (costStart.tv_sec * 1000000 + costStart.tv_usec)) / 1000000.0;
  totalCost += elapsed;
}

int main(int argc, char *argv[])
{
  char oriFilePath[640], outputFilePath[645];
  if (argc < 4)
  {
    printf("Usage: testfloat_compress_fastmode1 [srcFilePath] [block size] [err bound] [data_type]\n");
    printf("Example: testfloat_compress_fastmode1 testfloat_8_8_128.dat 64 1E-3 float\n");
    exit(0);
  }

  sprintf(oriFilePath, "%s", argv[1]);
  int blockSize = atoi(argv[2]);
  double errBound = atof(argv[3]);
  char* dataType = argv[4];
  int sz_dataType = SZ_FLOAT;
  if (strcmp(dataType, "float") != 0 && strcmp(dataType, "double") != 0)
  {
    printf("Error: data type must be 'float' or 'double'!\n");
    exit(0);
  }
  if (strcmp(dataType, "double") == 0)
  {
    sz_dataType = SZ_DOUBLE;
  }
  
  sprintf(outputFilePath, "%s.szp", oriFilePath);

  int status = 0;
  size_t nbEle;
  void *data;
  if (sz_dataType == SZ_FLOAT)
  {
    data = szp_readFloatData(oriFilePath, &nbEle, &status);
  }
  else // SZ_DOUBLE
  {
    data = szp_readDoubleData(oriFilePath, &nbEle, &status);
  }
  if (status != SZ_SCES)
  {
    printf("Error: data file %s cannot be read!\n", oriFilePath);
    exit(0);
  }


  size_t outSize;
  cost_start();
  unsigned char *bytes = szp_compress(SZP_NONRANDOMACCESS, sz_dataType, data, &outSize, ABS, errBound, 0, nbEle, blockSize);
  // unsigned char *bytes = szp_compress(SZP_RANDOMACCESS, SZ_FLOAT, data, &outSize, ABS, errBound, 0, nbEle, blockSize);
  cost_end();
  printf("\ntimecost=%f, total fastmode1\n", totalCost);
  printf("compression size = %zu, CR = %f, writing to %s\n", outSize, 1.0f * nbEle * ((sz_dataType == SZ_FLOAT) ? sizeof(float) : sizeof(double))/ outSize, outputFilePath);
  szp_writeByteData(bytes, outSize, outputFilePath, &status);
  if (status != SZ_SCES)
  {
    printf("Error: data file %s cannot be written!\n", outputFilePath);
    exit(0);
  }

  printf("done\n");
  free(bytes);
  free(data);

  return 0;
}
