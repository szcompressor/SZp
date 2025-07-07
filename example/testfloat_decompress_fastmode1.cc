/**
 *  @file testfloat_decompress_fastmode1.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using Decompression interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "szp.h"
#include <sys/time.h>
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
    size_t nbEle, totalNbEle;
    char zipFilePath[640], outputFilePath[645];
    if (argc < 5) // Corrected argument count check
    {
        printf("Usage: testfloat_decompress_fastmode1 [srcFilePath] [nbEle] [block size] [data_type]\n");
        printf("Example: testfloat_decompress_fastmode1 testfloat_8_8_128.dat.szp 8192 64 float\n");
        exit(0);
    }

    sprintf(zipFilePath, "%s", argv[1]);
    nbEle = strtoimax(argv[2], NULL, 10);
    int blockSize = atoi(argv[3]);
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

    sprintf(outputFilePath, "%s.out", zipFilePath);

    size_t byteLength;
    int status;
    unsigned char *bytes = szp_readByteData(zipFilePath, &byteLength, &status);
    if (status != SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", zipFilePath);
        exit(0);
    }

    cost_start();

    void *data = szp_decompress(SZP_RANDOMACCESS, sz_dataType, bytes, byteLength, nbEle, blockSize);
    // void *data = szp_decompress(SZP_NONRANDOMACCESS, sz_dataType, bytes, byteLength, nbEle, blockSize);
    cost_end();
    
    free(bytes);
    printf("timecost=%f\n", totalCost);

    // Write the decompressed data to the output file
    if (sz_dataType == SZ_FLOAT)
    {
        szp_writeFloatData_inBytes((float*)data, nbEle, outputFilePath, &status);
    }
    else // SZ_DOUBLE
    {
        szp_writeDoubleData_inBytes((double*)data, nbEle, outputFilePath, &status);
    }

    if (status != SZ_SCES)
    {
        printf("Error: %s cannot be written!\n", outputFilePath);
        exit(0);
    }
    printf("done\n");

    char oriFilePath[645];
    strcpy(oriFilePath, zipFilePath);
    oriFilePath[strlen(zipFilePath) - 4] = '\0';
    
    void* ori_data = NULL;
    if (sz_dataType == SZ_FLOAT)
    {
        ori_data = szp_readFloatData(oriFilePath, &totalNbEle, &status);
    }
    else // SZ_DOUBLE
    {
        ori_data = szp_readDoubleData(oriFilePath, &totalNbEle, &status);
    }
    
    if (status != SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", oriFilePath);
        free(data); // Free allocated memory before exiting
        exit(0);
    }

    if (sz_dataType == SZ_FLOAT)
    {
        float* ori_data_f = (float*)ori_data;
        float* data_f = (float*)data;
        size_t i = 0;
        double Max = 0, Min = 0, diffMax = 0; // Use double for precision in calculations

        if (nbEle > 0) {
            Max = ori_data_f[0];
            Min = ori_data_f[0];
            diffMax = fabs(data_f[0] - ori_data_f[0]);
        }

        double sum1 = 0, sum2 = 0;
        for (i = 0; i < nbEle; i++)
        {
            sum1 += ori_data_f[i];
            sum2 += data_f[i];
        }
        double mean1 = (nbEle == 0) ? 0 : sum1 / nbEle;
        double mean2 = (nbEle == 0) ? 0 : sum2 / nbEle;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0;
        double maxpw_relerr = 0;
        for (i = 0; i < nbEle; i++)
        {
            if (Max < ori_data_f[i]) Max = ori_data_f[i];
            if (Min > ori_data_f[i]) Min = ori_data_f[i];

            double err = fabs(data_f[i] - ori_data_f[i]);
            if (ori_data_f[i] != 0)
            {
                double relerr = err / fabs(ori_data_f[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }
            if (diffMax < err) diffMax = err;

            prodSum += (ori_data_f[i] - mean1) * (data_f[i] - mean2);
            sum3 += (ori_data_f[i] - mean1) * (ori_data_f[i] - mean1);
            sum4 += (data_f[i] - mean2) * (data_f[i] - mean2);
            sum += err * err;
        }

        double std1 = (nbEle == 0) ? 0 : sqrt(sum3 / nbEle);
        double std2 = (nbEle == 0) ? 0 : sqrt(sum4 / nbEle);
        double ee = (nbEle == 0) ? 0 : prodSum / nbEle;
        double acEff = (std1 * std2 == 0) ? 0 : ee / (std1 * std2);
        double mse = (nbEle == 0) ? 0 : sum / nbEle;
        double range = Max - Min;
        double psnr = (mse == 0) ? 999 : 20 * log10(range) - 10 * log10(mse);
        double nrmse = (range == 0) ? 0 : sqrt(mse) / range;
        double compressionRatio = 1.0 * nbEle * sizeof(float) / byteLength;

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.10f\n", diffMax);
        printf("Max relative error = %f\n", (range == 0) ? 0 : diffMax / range);
        printf("Max pw relative error = %f\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.20G\n", psnr, nrmse);
        printf("acEff=%f\n", acEff);
        printf("compressionRatio = %f\n", compressionRatio);
    }
    else // SZ_DOUBLE
    {
        double* ori_data_d = (double*)ori_data;
        double* data_d = (double*)data;
        size_t i = 0;
        double Max = 0, Min = 0, diffMax = 0;
        
        if (nbEle > 0) {
            Max = ori_data_d[0];
            Min = ori_data_d[0];
            diffMax = fabs(data_d[0] - ori_data_d[0]);
        }

        double sum1 = 0, sum2 = 0;
        for (i = 0; i < nbEle; i++)
        {
            sum1 += ori_data_d[i];
            sum2 += data_d[i];
        }
        double mean1 = (nbEle == 0) ? 0 : sum1 / nbEle;
        double mean2 = (nbEle == 0) ? 0 : sum2 / nbEle;

        double sum3 = 0, sum4 = 0;
        double sum = 0, prodSum = 0;
        double maxpw_relerr = 0;
        for (i = 0; i < nbEle; i++)
        {
            if (Max < ori_data_d[i]) Max = ori_data_d[i];
            if (Min > ori_data_d[i]) Min = ori_data_d[i];
            
            double err = fabs(data_d[i] - ori_data_d[i]);
            if (ori_data_d[i] != 0)
            {
                double relerr = err / fabs(ori_data_d[i]);
                if (maxpw_relerr < relerr)
                    maxpw_relerr = relerr;
            }
            if (diffMax < err) diffMax = err;

            prodSum += (ori_data_d[i] - mean1) * (data_d[i] - mean2);
            sum3 += (ori_data_d[i] - mean1) * (ori_data_d[i] - mean1);
            sum4 += (data_d[i] - mean2) * (data_d[i] - mean2);
            sum += err * err;
        }

        double std1 = (nbEle == 0) ? 0 : sqrt(sum3 / nbEle);
        double std2 = (nbEle == 0) ? 0 : sqrt(sum4 / nbEle);
        double ee = (nbEle == 0) ? 0 : prodSum / nbEle;
        double acEff = (std1 * std2 == 0) ? 0 : ee / (std1 * std2);
        double mse = (nbEle == 0) ? 0 : sum / nbEle;
        double range = Max - Min;
        double psnr = (mse == 0) ? 999 : 20 * log10(range) - 10 * log10(mse);
        double nrmse = (range == 0) ? 0 : sqrt(mse) / range;
        double compressionRatio = 1.0 * nbEle * sizeof(double) / byteLength;

        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.10f\n", diffMax);
        printf("Max relative error = %f\n", (range == 0) ? 0 : diffMax / range);
        printf("Max pw relative error = %f\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.20G\n", psnr, nrmse);
        printf("acEff=%f\n", acEff);
        printf("compressionRatio = %f\n", compressionRatio);
    }

    free(data);
    free(ori_data);
    return 0;
}
