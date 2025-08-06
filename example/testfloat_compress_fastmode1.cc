/**
 *  @file testfloat_compress_fastmode1.cc
 *  @author Jiajun Huang, Sheng Di (modified by Jeffey Chou)
 *  @date Oct, 2023 (modified Aug, 2025)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "szp.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include <sys/time.h>
#include <vector>
#include <algorithm>

struct timeval startTime;
struct timeval endTime;   /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;

// Structure to hold benchmark results
typedef struct {
    double timeCost;
    size_t compressedSize;
    double compressionRatio;
} BenchmarkResult;

// Function to write results to CSV file
void writeToCsv(const char* filename, const std::vector<BenchmarkResult>& results) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        printf("Error: Cannot open file %s for writing\n", filename);
        return;
    }
    
    fprintf(fp, "Run,TimeSeconds,CompressionSize,CompressionRatio\n");
    for (size_t i = 0; i < results.size(); i++) {
        fprintf(fp, "%zu,%f,%zu,%f\n", i+1, results[i].timeCost, 
                results[i].compressedSize, results[i].compressionRatio);
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
    
    // CR statistics
    double crSum = 0, crMin = results[0].compressionRatio, crMax = results[0].compressionRatio;
    for (const auto& res : results) {
        crSum += res.compressionRatio;
        crMin = std::min(crMin, res.compressionRatio);
        crMax = std::max(crMax, res.compressionRatio);
    }
    double crAvg = crSum / results.size();
    
    double crVariance = 0;
    for (const auto& res : results) {
        crVariance += (res.compressionRatio - crAvg) * (res.compressionRatio - crAvg);
    }
    double crStdDev = sqrt(crVariance / results.size());
    
    // Print statistics
    printf("\n======= Performance Statistics (excluding warmup runs) =======\n");
    printf("Time (seconds): Min=%.6f, Max=%.6f, Avg=%.6f, StdDev=%.6f\n", 
           timeMin, timeMax, timeAvg, timeStdDev);
    printf("Compression Ratio: Min=%.6f, Max=%.6f, Avg=%.6f, StdDev=%.6f\n", 
           crMin, crMax, crAvg, crStdDev);
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

void printUsage() {
    printf("Usage: testfloat_compress_fastmode1 [srcFilePath] [block size] [err bound] [data_type] [access_type] [repetitions] [warmup_runs] [output_csv]\n");
    printf("Example: testfloat_compress_fastmode1 testfloat_8_8_128.dat 64 1E-3 float random 10 2 results.csv\n");
    printf("  - repetitions, warmup_runs, and output_csv are optional\n");
}

int main(int argc, char *argv[])
{
    char oriFilePath[640], outputFilePath[645], csvFilePath[650] = {0};
    if (argc < 6)
    {
        printUsage();
        exit(0);
    }

    sprintf(oriFilePath, "%s", argv[1]);
    int blockSize = atoi(argv[2]);
    double errBound = atof(argv[3]);
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

    sprintf(outputFilePath, "%s.szp", oriFilePath);

    // Print benchmark configuration
    printf("SZp Compression Benchmark Configuration\n");
    printf("======================================\n");
    printf("Input file: %s\n", oriFilePath);
    printf("Block size: %d\n", blockSize);
    printf("Error bound: %g\n", errBound);
    printf("Data type: %s\n", dataType);
    printf("Access type: %s\n", accessType);
    printf("Repetitions: %d\n", repetitions);
    printf("Warmup runs: %d\n", warmupRuns);
    if (csvFilePath[0]) printf("CSV output: %s\n", csvFilePath);
    printf("======================================\n\n");

    // Read data only once
    int status = 0;
    size_t nbEle;
    void *data;
    printf("Reading input data...\n");
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
    printf("Data elements: %zu\n", nbEle);

    // Vector to store benchmark results (excluding warmup runs)
    std::vector<BenchmarkResult> results;

    // Perform warmup runs (results not recorded)
    if (warmupRuns > 0) {
        printf("\nPerforming %d warmup runs...\n", warmupRuns);
        for (int i = 0; i < warmupRuns; i++) {
            printf("  Warmup run %d of %d: ", i+1, warmupRuns);
            
            size_t outSize;
            cost_start();
            unsigned char *bytes = szp_compress(accessMode, sz_dataType, data, &outSize, 
                                               ABS, errBound, 0, nbEle, blockSize);
            cost_end();
            
            printf("time=%f, size=%zu, CR=%f\n", 
                   totalCost, outSize, 
                   1.0f * nbEle * ((sz_dataType == SZ_FLOAT) ? sizeof(float) : sizeof(double)) / outSize);
            
            free(bytes);
        }
    }

    // Perform actual benchmark runs
    printf("\nStarting %d benchmark runs...\n", repetitions);
    for (int i = 0; i < repetitions; i++) {
        printf("\nRun %d of %d:\n", i+1, repetitions);
        
        size_t outSize;
        cost_start();
        unsigned char *bytes = szp_compress(accessMode, sz_dataType, data, &outSize, 
                                           ABS, errBound, 0, nbEle, blockSize);
        cost_end();
        
        double compressionRatio = 1.0f * nbEle * ((sz_dataType == SZ_FLOAT) ? sizeof(float) : sizeof(double)) / outSize;
        printf("  timecost=%f, compression size=%zu, CR=%f\n", totalCost, outSize, compressionRatio);
        
        // Save result
        BenchmarkResult result;
        result.timeCost = totalCost;
        result.compressedSize = outSize;
        result.compressionRatio = compressionRatio;
        results.push_back(result);
        
        // Write compressed data to file only on the last run
        if (i == repetitions - 1) {
            printf("  Writing compressed data to %s\n", outputFilePath);
            szp_writeByteData(bytes, outSize, outputFilePath, &status);
            if (status != SZ_SCES) {
                printf("Error: data file %s cannot be written!\n", outputFilePath);
            }
        }
        
        free(bytes);
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
    free(data);

    return 0;
}