/**
 *  @file test_multidim.cc
 *  @author Jiefeng Zhou
 *  @date Aug, 2025
 *  @brief A benchmark tool for compressing and decompressing multi-dimensional datasets.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "szp.h"
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <numeric>

// --- Utility Functions ---

double totalCost = 0;
struct timeval costStart;

void cost_start() {
    gettimeofday(&costStart, NULL);
}

void cost_end() {
    struct timeval costEnd;
    gettimeofday(&costEnd, NULL);
    totalCost = ((costEnd.tv_sec * 1000000 + costEnd.tv_usec) - (costStart.tv_sec * 1000000 + costStart.tv_usec)) / 1000000.0;
}

// --- Result Structures and Analysis ---

struct CompResult {
    double timeCost;
    size_t compressedSize;
    double compressionRatio;
};

struct DecompResult {
    double timeCost;
    double maxAbsErr;
    double psnr;
    double nrmse;
    double compressionRatio;
};

// --- Compression Statistics ---
void writeCompToCsv(const char* filename, const std::vector<CompResult>& results) {
    FILE* fp = fopen(filename, "w");
    if (!fp) return;
    fprintf(fp, "Run,TimeSeconds,CompressedSize,CompressionRatio\n");
    for (size_t i = 0; i < results.size(); i++) {
        fprintf(fp, "%zu,%.6f,%zu,%.2f\n", i + 1, results[i].timeCost, results[i].compressedSize, results[i].compressionRatio);
    }
    fclose(fp);
}

void calculateCompStats(const std::vector<CompResult>& results) {
    if (results.empty()) return;
    std::vector<double> times, crs;
    for(const auto& res : results) {
        times.push_back(res.timeCost);
        crs.push_back(res.compressionRatio);
    }
    
    double timeAvg = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double crAvg = std::accumulate(crs.begin(), crs.end(), 0.0) / crs.size();
    
    double timeStdDev = 0.0;
    for(const auto& t : times) timeStdDev += (t - timeAvg) * (t - timeAvg);
    timeStdDev = sqrt(timeStdDev / times.size());

    double crStdDev = 0.0;
    for(const auto& cr : crs) crStdDev += (cr - crAvg) * (cr - crAvg);
    crStdDev = sqrt(crStdDev / crs.size());

    printf("\n======= Compression Performance (over %zu runs) =======\n", results.size());
    printf("Time (s):   Avg=%.6f, StdDev=%.6f\n", timeAvg, timeStdDev);
    printf("Comp Ratio: Avg=%.2f, StdDev=%.2f\n", crAvg, crStdDev);
    printf("======================================================\n");
}


// --- Decompression Statistics ---
void writeDecompToCsv(const char* filename, const std::vector<DecompResult>& results) {
    FILE* fp = fopen(filename, "w");
    if (!fp) return;
    fprintf(fp, "Run,TimeSeconds,MaxAbsErr,PSNR,NRMSE,CompressionRatio\n");
    for (size_t i = 0; i < results.size(); i++) {
        fprintf(fp, "%zu,%.6f,%.6e,%.4f,%.6e,%.2f\n", i + 1, results[i].timeCost, results[i].maxAbsErr, results[i].psnr, results[i].nrmse, results[i].compressionRatio);
    }
    fclose(fp);
}

void calculateDecompStats(const std::vector<DecompResult>& results) {
    if (results.empty()) return;
    std::vector<double> times, psnrs;
    for(const auto& res : results) {
        times.push_back(res.timeCost);
        psnrs.push_back(res.psnr);
    }

    double timeAvg = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double psnrAvg = std::accumulate(psnrs.begin(), psnrs.end(), 0.0) / psnrs.size();

    double timeStdDev = 0.0;
    for(const auto& t : times) timeStdDev += (t - timeAvg) * (t - timeAvg);
    timeStdDev = sqrt(timeStdDev / times.size());

    double psnrStdDev = 0.0;
    for(const auto& p : psnrs) psnrStdDev += (p - psnrAvg) * (p - psnrAvg);
    psnrStdDev = sqrt(psnrStdDev / psnrs.size());

    printf("\n======= Decompression Performance (over %zu runs) =======\n", results.size());
    printf("Time (s): Avg=%.6f, StdDev=%.6f\n", timeAvg, timeStdDev);
    printf("PSNR:     Avg=%.4f, StdDev=%.4f\n", psnrAvg, psnrStdDev);
    printf("=========================================================\n");
}


template<typename T>
DecompResult evaluate(const T* ori_data, const T* dec_data, size_t num_elements, size_t compressed_size) {
    DecompResult res;
    double max_val = ori_data[0], min_val = ori_data[0], max_abs_err = 0.0, mse = 0.0;

    for (size_t i = 0; i < num_elements; ++i) {
        if (ori_data[i] > max_val) max_val = ori_data[i];
        if (ori_data[i] < min_val) min_val = ori_data[i];
        double abs_err = fabs((double)ori_data[i] - (double)dec_data[i]);
        if (abs_err > max_abs_err) max_abs_err = abs_err;
        mse += abs_err * abs_err;
    }
    mse /= num_elements;
    double value_range = max_val - min_val;
    
    res.maxAbsErr = max_abs_err;
    res.psnr = (mse == 0) ? 999.99 : 20 * log10(value_range) - 10 * log10(mse);
    res.nrmse = (value_range == 0) ? 0 : sqrt(mse) / value_range;
    res.compressionRatio = (double)(num_elements * sizeof(T)) / compressed_size;

    printf("  Max Abs Error: %.6e, PSNR: %.4f, NRMSE: %.6e, CR: %.2f\n", res.maxAbsErr, res.psnr, res.nrmse, res.compressionRatio);
    return res;
}

// --- Main Logic ---

void print_usage() {
    printf("Usage: test_multidim <mode> <options>\n");
    printf("\n  Mode: -c (compress) or -d (decompress)\n");
    printf("\n  Compress Mode:\n");
    printf("    ./test_multidim -c <filepath> <dtype> <err_bound> <block_size> <dims...> [-r reps] [-w warmups] [-o out.csv]\n");
    printf("    Example: ./test_multidim -c temp.f32 float 1e-4 128 512 512 512 -r 10 -w 2\n");
    printf("\n  Decompress Mode:\n");
    printf("    ./test_multidim -d <compressed_file> <dtype> <block_size> <dims...> [-r reps] [-w warmups] [-o out.csv]\n");
    printf("    Example: ./test_multidim -d temp.f32.szp float 128 512 512 512\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    bool compress_mode = (strcmp(argv[1], "-c") == 0);
    bool decompress_mode = (strcmp(argv[1], "-d") == 0);

    if (!compress_mode && !decompress_mode) {
        print_usage();
        return 1;
    }

    // --- Argument Parsing ---
    const char* filepath = argv[2];
    const char* datatype_str = argv[3];
    double err_bound = compress_mode ? atof(argv[4]) : 0;
    int block_size = compress_mode ? atoi(argv[5]) : atoi(argv[4]);
    
    int dim_start_index = compress_mode ? 6 : 5;
    size_t num_elements = 1;
    int current_arg = dim_start_index;
    while (current_arg < argc && argv[current_arg][0] != '-') {
        num_elements *= (size_t)atol(argv[current_arg]);
        current_arg++;
    }

    int repetitions = 1, warmupRuns = 0;
    char csvFilePath[256] = {0};
    while (current_arg < argc) {
        if (strcmp(argv[current_arg], "-r") == 0) repetitions = atoi(argv[++current_arg]);
        else if (strcmp(argv[current_arg], "-w") == 0) warmupRuns = atoi(argv[++current_arg]);
        else if (strcmp(argv[current_arg], "-o") == 0) sprintf(csvFilePath, "%s", argv[++current_arg]);
        current_arg++;
    }

    int sz_datatype = (strcmp(datatype_str, "float") == 0) ? SZ_FLOAT : SZ_DOUBLE;
    size_t dtype_size = (sz_datatype == SZ_FLOAT) ? sizeof(float) : sizeof(double);

    // --- Compression Logic ---
    if (compress_mode) {
        int status;
        void* data = (sz_datatype == SZ_FLOAT) ? (void*)szp_readFloatData((char*)filepath, &num_elements, &status) : (void*)szp_readDoubleData((char*)filepath, &num_elements, &status);
        if (status != SZ_SCES) {
            fprintf(stderr, "Error reading %s\n", filepath);
            return 1;
        }

        printf("--- Compressing %s ---\n", filepath);
        std::vector<CompResult> results;
        for (int i = 0; i < warmupRuns + repetitions; ++i) {
            size_t compressed_size;
            cost_start();
            unsigned char* bytes = szp_compress(SZP_RANDOMACCESS, sz_datatype, data, &compressed_size, ABS, err_bound, 0, num_elements, block_size);
            cost_end();

            if (i < warmupRuns) {
                printf("Warmup %d/%d: Time=%.4fs\n", i + 1, warmupRuns, totalCost);
            } else {
                printf("Run %d/%d: Time=%.4fs, Size=%zu, CR=%.2f\n", i - warmupRuns + 1, repetitions, totalCost, compressed_size, (double)(num_elements * dtype_size) / compressed_size);
                CompResult res = {totalCost, compressed_size, (double)(num_elements * dtype_size) / compressed_size};
                results.push_back(res);
                if (i == warmupRuns + repetitions - 1) {
                    char outpath[256];
                    sprintf(outpath, "%s.szp", filepath);
                    szp_writeByteData(bytes, compressed_size, outpath, &status);
                    printf("Compressed data written to %s\n", outpath);
                }
            }
            free(bytes);
        }
        free(data);
        if (results.size() > 1) calculateCompStats(results);
        if (csvFilePath[0]) writeCompToCsv(csvFilePath, results);
    }

    // --- Decompression Logic ---
    if (decompress_mode) {
        int status;
        size_t compressed_size;
        unsigned char* bytes = szp_readByteData((char*)filepath, &compressed_size, &status);
        if (status != SZ_SCES) {
            fprintf(stderr, "Error reading %s\n", filepath);
            return 1;
        }

        char oriFilePath[256];
        strcpy(oriFilePath, filepath);
        oriFilePath[strlen(filepath) - 4] = '\0'; // Remove .szp
        void* ori_data = (sz_datatype == SZ_FLOAT) ? (void*)szp_readFloatData(oriFilePath, &num_elements, &status) : (void*)szp_readDoubleData(oriFilePath, &num_elements, &status);
        if (status != SZ_SCES) {
            fprintf(stderr, "Error reading original file %s for verification\n", oriFilePath);
            free(bytes);
            return 1;
        }

        printf("--- Decompressing %s ---\n", filepath);
        std::vector<DecompResult> results;
        for (int i = 0; i < warmupRuns + repetitions; ++i) {
            cost_start();
            void* dec_data = szp_decompress(SZP_RANDOMACCESS, sz_datatype, bytes, compressed_size, num_elements, block_size);
            cost_end();

            if (i < warmupRuns) {
                printf("Warmup %d/%d: Time=%.4fs\n", i + 1, warmupRuns, totalCost);
            } else {
                printf("Run %d/%d: Time=%.4fs\n", i - warmupRuns + 1, repetitions, totalCost);
                DecompResult res = (sz_datatype == SZ_FLOAT) ? evaluate((float*)ori_data, (float*)dec_data, num_elements, compressed_size) : evaluate((double*)ori_data, (double*)dec_data, num_elements, compressed_size);
                res.timeCost = totalCost;
                results.push_back(res);
                if (i == warmupRuns + repetitions - 1) {
                    char outpath[256];
                    sprintf(outpath, "%s.out", filepath);
                    if(sz_datatype == SZ_FLOAT) szp_writeFloatData_inBytes((float*)dec_data, num_elements, outpath, &status);
                    else szp_writeDoubleData_inBytes((double*)dec_data, num_elements, outpath, &status);
                    printf("Decompressed data written to %s\n", outpath);
                }
            }
            free(dec_data);
        }
        free(bytes);
        free(ori_data);
        if (results.size() > 1) calculateDecompStats(results);
        if (csvFilePath[0]) writeDecompToCsv(csvFilePath, results);
    }

    return 0;
}