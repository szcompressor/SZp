/**
 *  @file testfloat_decompress_fastmode1.cc
 *  @author Sheng Di (modified by Jeffey Chou)
 *  @date April, 2015 (modified Aug, 2025)
 *  @brief This is an example of using Decompression interface with integrated benchmarking.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "szp.h"
#include <sys/time.h>
#ifdef _OPENMP
#include "omp.h"
#endif
#include <vector>
#include <algorithm>

struct timeval startTime;
struct timeval endTime;   /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

// Structure to hold benchmark results
typedef struct {
    double timeCost;
    double maxAbsErr;
    double maxRelErr;
    double maxPwRelErr;
    double psnr;
    double nrmse;
    double acEff;
    double compressionRatio;
} BenchmarkResult;

// Function to write results to CSV file
void writeToCsv(const char* filename, const std::vector<BenchmarkResult>& results) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        printf("Error: Cannot open file %s for writing\n", filename);
        return;
    }
    
    fprintf(fp, "Run,TimeSeconds,MaxAbsErr,MaxRelErr,MaxPwRelErr,PSNR,NRMSE,AccEff,CompressionRatio\n");
    for (size_t i = 0; i < results.size(); i++) {
        fprintf(fp, "%zu,%f,%f,%f,%f,%f,%e,%f,%f\n", i+1, 
                results[i].timeCost, results[i].maxAbsErr, results[i].maxRelErr,
                results[i].maxPwRelErr, results[i].psnr, results[i].nrmse,
                results[i].acEff, results[i].compressionRatio);
    }
    fclose(fp);
}

// Function to calculate statistics on benchmark results
void calculateStats(const std::vector<BenchmarkResult>& results) {
    if (results.empty()) return;
    
    // Time statistics
    double timeSum = 0, timeMin = results[0].timeCost, timeMax = results[0].timeCost;
    for (const auto& res : results) {
        timeSum += res.timeCost;
        timeMin = std::min(timeMin, res.timeCost);
        timeMax = std::max(timeMax, res.timeCost);
    }
    double timeAvg = timeSum / results.size();
    
    double timeVariance = 0;
    for (const auto& res : results) {
        timeVariance += (res.timeCost - timeAvg) * (res.timeCost - timeAvg);
    }
    double timeStdDev = sqrt(timeVariance / results.size());
    
    // PSNR statistics
    double psnrSum = 0, psnrMin = results[0].psnr, psnrMax = results[0].psnr;
    for (const auto& res : results) {
        psnrSum += res.psnr;
        psnrMin = std::min(psnrMin, res.psnr);
        psnrMax = std::max(psnrMax, res.psnr);
    }
    double psnrAvg = psnrSum / results.size();
    
    double psnrVariance = 0;
    for (const auto& res : results) {
        psnrVariance += (res.psnr - psnrAvg) * (res.psnr - psnrAvg);
    }
    double psnrStdDev = sqrt(psnrVariance / results.size());
    
    // Print statistics
    printf("\n======= Performance Statistics (excluding warmup runs) =======\n");
    printf("Time (seconds): Min=%.6f, Max=%.6f, Avg=%.6f, StdDev=%.6f\n", 
           timeMin, timeMax, timeAvg, timeStdDev);
    printf("PSNR: Min=%.6f, Max=%.6f, Avg=%.6f, StdDev=%.6f\n", 
           psnrMin, psnrMax, psnrAvg, psnrStdDev);
    printf("===========================================================\n");
}

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

// Function to evaluate decompression quality
BenchmarkResult evaluateQuality(void* ori_data, void* data, size_t nbEle, int sz_dataType, size_t byteLength) {
    BenchmarkResult result;
    result.timeCost = totalCost;
    
    if (sz_dataType == SZ_FLOAT)
    {
        float* ori_data_f = (float*)ori_data;
        float* data_f = (float*)data;
        size_t i = 0;
        double Max = 0, Min = 0, diffMax = 0;

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

        result.maxAbsErr = diffMax;
        result.maxRelErr = (range == 0) ? 0 : diffMax / range;
        result.maxPwRelErr = maxpw_relerr;
        result.psnr = psnr;
        result.nrmse = nrmse;
        result.acEff = acEff;
        result.compressionRatio = compressionRatio;
        
        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.10f\n", diffMax);
        printf("Max relative error = %f\n", result.maxRelErr);
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

        result.maxAbsErr = diffMax;
        result.maxRelErr = (range == 0) ? 0 : diffMax / range;
        result.maxPwRelErr = maxpw_relerr;
        result.psnr = psnr;
        result.nrmse = nrmse;
        result.acEff = acEff;
        result.compressionRatio = compressionRatio;
        
        printf("Min=%.20G, Max=%.20G, range=%.20G\n", Min, Max, range);
        printf("Max absolute error = %.10f\n", diffMax);
        printf("Max relative error = %f\n", result.maxRelErr);
        printf("Max pw relative error = %f\n", maxpw_relerr);
        printf("PSNR = %f, NRMSE= %.20G\n", psnr, nrmse);
        printf("acEff=%f\n", acEff);
        printf("compressionRatio = %f\n", compressionRatio);
    }
    
    return result;
}

void printUsage() {
    printf("Usage: testfloat_decompress_fastmode1 [srcFilePath] [nbEle] [block size] [data_type] [access_type] [repetitions] [warmup_runs] [output_csv]\n");
    printf("Example: testfloat_decompress_fastmode1 testfloat_8_8_128.dat.szp 8192 64 float random 10 2 results.csv\n");
    printf("  - repetitions, warmup_runs, and output_csv are optional\n");
}

int main(int argc, char *argv[])
{
    size_t nbEle, totalNbEle;
    char zipFilePath[640], outputFilePath[645], csvFilePath[650] = {0};
    
    if (argc < 6)
    {
        printUsage();
        exit(0);
    }

    sprintf(zipFilePath, "%s", argv[1]);
    nbEle = strtoimax(argv[2], NULL, 10);
    int blockSize = atoi(argv[3]);
    char* dataType = argv[4];
    char* accessType = argv[5];
    
    // Benchmark parameters - default values
    int repetitions = 1;
    int warmupRuns = 0;
    
    // Parse optional benchmark parameters
    if (argc >= 7) repetitions = atoi(argv[6]);
    if (argc >= 8) warmupRuns = atoi(argv[7]);
    if (argc >= 9) sprintf(csvFilePath, "%s", argv[8]);
    
    int sz_dataType = SZ_FLOAT;
    int accessMode = SZP_RANDOMACCESS;

    if (strcmp(dataType, "float") != 0 && strcmp(dataType, "double") != 0)
    {
        printf("Error: data type must be 'float' or 'double'!\n");
        exit(0);
    }
    if (strcmp(dataType, "double") == 0)
    {
        sz_dataType = SZ_DOUBLE;
    }

    if (strcmp(accessType, "random") == 0)
    {
        accessMode = SZP_RANDOMACCESS;
    }
    else if (strcmp(accessType, "nonrandom") == 0)
    {
        accessMode = SZP_NONRANDOMACCESS;
    }
    else
    {
        printf("Error: access type must be 'random' or 'nonrandom'!\n");
        exit(0);
    }

    // Print benchmark configuration
    printf("SZp Decompression Benchmark Configuration\n");
    printf("=========================================\n");
    printf("Compressed file: %s\n", zipFilePath);
    printf("Data elements: %zu\n", nbEle);
    printf("Block size: %d\n", blockSize);
    printf("Data type: %s\n", dataType);
    printf("Access type: %s\n", accessType);
    printf("Repetitions: %d\n", repetitions);
    printf("Warmup runs: %d\n", warmupRuns);
    if (csvFilePath[0]) printf("CSV output: %s\n", csvFilePath);
    printf("=========================================\n\n");

    sprintf(outputFilePath, "%s.out", zipFilePath);

    // Read compressed data only once
    size_t byteLength;
    int status;
    printf("Reading compressed file...\n");
    unsigned char *bytes = szp_readByteData(zipFilePath, &byteLength, &status);
    if (status != SZ_SCES)
    {
        printf("Error: %s cannot be read!\n", zipFilePath);
        exit(0);
    }

    // Read original data for quality evaluation
    char oriFilePath[645];
    strcpy(oriFilePath, zipFilePath);
    oriFilePath[strlen(zipFilePath) - 4] = '\0';
    
    printf("Reading original file for comparison...\n");
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
        free(bytes);
        exit(0);
    }

    if (nbEle > totalNbEle) {
        printf("Warning: Specified element count (%zu) exceeds file size (%zu). Using file size.\n", 
               nbEle, totalNbEle);
        nbEle = totalNbEle;
    }

    // Vector to store benchmark results (excluding warmup runs)
    std::vector<BenchmarkResult> results;

    // Perform warmup runs (results not recorded)
    if (warmupRuns > 0) {
        printf("\nPerforming %d warmup runs...\n", warmupRuns);
        for (int i = 0; i < warmupRuns; i++) {
            printf("  Warmup run %d of %d: ", i+1, warmupRuns);
            
            cost_start();
            void *data = szp_decompress(accessMode, sz_dataType, bytes, byteLength, nbEle, blockSize);
            cost_end();
            
            printf("time=%f\n", totalCost);
            free(data);
        }
    }

    // Perform actual benchmark runs
    printf("\nStarting %d benchmark runs...\n", repetitions);
    for (int i = 0; i < repetitions; i++) {
        printf("\nRun %d of %d:\n", i+1, repetitions);
        
        cost_start();
        void *data = szp_decompress(accessMode, sz_dataType, bytes, byteLength, nbEle, blockSize);
        cost_end();
        
        printf("timecost=%f\n", totalCost);
        
        // Evaluate decompression quality
        BenchmarkResult result = evaluateQuality(ori_data, data, nbEle, sz_dataType, byteLength);
        result.timeCost = totalCost;  // Set the measured time
        
        // Save result
        results.push_back(result);
        
        // Write decompressed data to file only on the last run
        if (i == repetitions - 1) {
            printf("  Writing decompressed data to %s\n", outputFilePath);
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
            }
        }
        
        free(data);
    }

    // Calculate and print statistics
    if (results.size() > 1) {
        calculateStats(results);
    }
    
    // Write to CSV if requested
    if (csvFilePath[0]) {
        writeToCsv(csvFilePath, results);
        printf("Results written to CSV: %s\n", csvFilePath);
    }

    printf("\nBenchmark completed.\n");
    free(bytes);
    free(ori_data);

    return 0;
}