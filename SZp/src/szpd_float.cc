/**
 *  @file szpd_float.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>, Sheng Di <sdi1@anl.gov>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "szpd_float.h"
#include <assert.h>
#include <math.h>
#include "szp_TypeManager.h"
#include "szp_CompressionToolkit.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace szp;

float *szp_float_decompress_openmp_threadblock(size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    float *newData = (float *)malloc(sizeof(float) * nbEle);
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    size_t threadblocksize = 0;
    int block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            rcp = cmpBytes + nbThreads * sizeof(size_t);
            threadblocksize = nbEle / nbThreads;
        }
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        float *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;

        int prior = 0;
        int current = 0;
        int diff = 0;

        unsigned int bit_count = 0;
        unsigned char *block_pointer = rcp + offsets[tid];

        float ori_prior = 0.0;
        float ori_current = 0.0;

        if (lo < hi) { // Ensure thread has data to process
            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);

            ori_prior = (float)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(float)); 
            newData_perthread += 1;
        }
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
        unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo + 1; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            bit_count = block_pointer[0];
            block_pointer++;
            
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                
                for (j = 0; j < current_block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                }
            }
            else
            {
                convertByteArray2IntArray_fast_1b_args(current_block_size, block_pointer, (current_block_size - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((current_block_size - 1) / 8 + 1);

                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, current_block_size, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;
                for (j = 0; j < current_block_size; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    newData_perthread++;
                }
            }
        }
        free(temp_predict_arr);
        free(temp_sign_arr);
    }
    return newData;

#else
    printf("Error! OpenMP not supported!\n");
    return NULL; 
#endif
}

void szp_float_decompress_single_thread_arg(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes)
{

    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    size_t threadblocksize = 0;
    
    int block_size = blockSize;
    int num_full_block_in_tb = 0;
    int num_remainder_in_tb = 0;

    nbThreads = 1;
    rcp = cmpBytes + nbThreads * sizeof(size_t);
    threadblocksize = nbEle / nbThreads;
    
    num_full_block_in_tb = (threadblocksize - 1) / block_size; 
    num_remainder_in_tb = (threadblocksize - 1) % block_size;

    int tid = 0;
    size_t lo = tid * threadblocksize;
    size_t hi = (tid + 1) * threadblocksize;
    float *newData_perthread = newData + lo;
    size_t i = 0;
    int j = 0;

    int prior = 0;
    int current = 0;
    int diff = 0;

    int bit_count = 0;
    unsigned char *outputBytes_perthread = rcp + offsets[tid]; 
    unsigned char *block_pointer = outputBytes_perthread;
    memcpy(&prior, block_pointer, sizeof(int));
    block_pointer += sizeof(unsigned int);

    float ori_prior = 0.0;
    float ori_current = 0.0;

    ori_prior = (float)prior * absErrBound;
    memcpy(newData_perthread, &ori_prior, sizeof(float)); 
    newData_perthread += 1;
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int savedbitsbytelength = 0;
    if (num_full_block_in_tb > 0)
    {
        for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
        {
            bit_count = block_pointer[0];
            block_pointer++;
            
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                
                for (j = 0; j < block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                }
            }
            else
            {
                
                convertByteArray2IntArray_fast_1b_args(block_size, block_pointer, (block_size - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((block_size - 1) / 8 + 1);

                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, block_size, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;
                for (j = 0; j < block_size; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    newData_perthread++;
                }
            }
        }
    }
    
    if (num_remainder_in_tb > 0)
    {
        for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
        {
            bit_count = block_pointer[0];
            block_pointer++;
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                }
            }
            else
            {
                convertByteArray2IntArray_fast_1b_args(num_remainder_in_tb, block_pointer, (num_remainder_in_tb - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((num_remainder_in_tb - 1) / 8 + 1);

                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, num_remainder_in_tb, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    newData_perthread++;
                }
            }
        }
    }
    free(temp_predict_arr);
    free(temp_sign_arr);
}

size_t szp_float_decompress_single_thread_arg_record(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes)
{
    size_t total_memaccess = 0;
    
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    size_t threadblocksize = 0;
    
    int block_size = blockSize;
    int num_full_block_in_tb = 0;
    int num_remainder_in_tb = 0;

    nbThreads = 1;
    rcp = cmpBytes + nbThreads * sizeof(size_t);
    threadblocksize = nbEle / nbThreads;
    
    num_full_block_in_tb = (threadblocksize - 1) / block_size; 
    num_remainder_in_tb = (threadblocksize - 1) % block_size;

    int tid = 0;
    size_t lo = tid * threadblocksize;
    size_t hi = (tid + 1) * threadblocksize;
    float *newData_perthread = newData + lo;
    size_t i = 0;
    int j = 0;

    int prior = 0;
    int current = 0;
    int diff = 0;

    int bit_count = 0;
    unsigned char *outputBytes_perthread = rcp + offsets[tid]; 
    total_memaccess += sizeof(size_t);
    unsigned char *block_pointer = outputBytes_perthread;
    memcpy(&prior, block_pointer, sizeof(int));
    block_pointer += sizeof(unsigned int);

    total_memaccess += sizeof(int);
    total_memaccess += sizeof(int);

    float ori_prior = 0.0;
    float ori_current = 0.0;

    ori_prior = (float)prior * absErrBound;
    memcpy(newData_perthread, &ori_prior, sizeof(float)); 
    total_memaccess += sizeof(float);
    total_memaccess += sizeof(float);
    newData_perthread += 1;
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int savedbitsbytelength = 0;
    if (num_full_block_in_tb > 0)
    {
        for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
        {
            bit_count = block_pointer[0];
            total_memaccess += sizeof(unsigned char);
            block_pointer++;
            
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                
                for (j = 0; j < block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                    total_memaccess += sizeof(float);
                    total_memaccess += sizeof(float);
                }
            }
            else
            {
                
                convertByteArray2IntArray_fast_1b_args(block_size, block_pointer, (block_size - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((block_size - 1) / 8 + 1);
                total_memaccess += (sizeof(unsigned int) * block_size);
                total_memaccess += (sizeof(unsigned char) * ((block_size - 1) / 8 + 1));
                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, block_size, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;
                total_memaccess += (sizeof(unsigned int) * block_size);
                total_memaccess += (sizeof(unsigned char) * savedbitsbytelength);
                for (j = 0; j < block_size; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                        total_memaccess += sizeof(unsigned int);
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                        total_memaccess += sizeof(unsigned int);
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    total_memaccess += sizeof(float);
                    total_memaccess += sizeof(float);
                    newData_perthread++;
                }
            }
        }
    }
    
    if (num_remainder_in_tb > 0)
    {
        for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
        {
            bit_count = block_pointer[0];
            total_memaccess += sizeof(unsigned char);
            block_pointer++;
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                    total_memaccess += sizeof(float);
                    total_memaccess += sizeof(float);
                }
            }
            else
            {
                convertByteArray2IntArray_fast_1b_args(num_remainder_in_tb, block_pointer, (num_remainder_in_tb - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((num_remainder_in_tb - 1) / 8 + 1);

                total_memaccess += (sizeof(unsigned int) * num_remainder_in_tb);
                total_memaccess += (sizeof(unsigned char) * ((num_remainder_in_tb - 1) / 8 + 1));

                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, num_remainder_in_tb, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;

                total_memaccess += (sizeof(unsigned int) * num_remainder_in_tb);
                total_memaccess += (sizeof(unsigned char) * savedbitsbytelength);
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                        total_memaccess += sizeof(unsigned int);
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                        total_memaccess += sizeof(unsigned int);
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    newData_perthread++;
                    total_memaccess += sizeof(float);
                    total_memaccess += sizeof(float);
                }
            }
        }
    }
    free(temp_predict_arr);
    free(temp_sign_arr);

    return total_memaccess;
}

void szp_float_decompress_openmp_threadblock_arg(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    size_t threadblocksize = 0;
    size_t block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            rcp = cmpBytes + nbThreads * sizeof(size_t);
            threadblocksize = nbEle / nbThreads;
        }
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        float *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;

        int prior = 0;
        int current = 0;
        int diff = 0;

        unsigned int bit_count = 0;
        unsigned char *block_pointer = rcp + offsets[tid];

        float ori_prior = 0.0;
        float ori_current = 0.0;

        if (lo < hi) { // Ensure thread has data to process
            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);

            ori_prior = (float)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(float)); 
            newData_perthread += 1;
        }
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
        unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo + 1; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            bit_count = block_pointer[0];
            block_pointer++;
            
            if (bit_count == 0)
            {
                ori_prior = (float)prior * absErrBound;
                
                for (j = 0; j < current_block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(float));
                    newData_perthread++;
                }
            }
            else
            {
                convertByteArray2IntArray_fast_1b_args(current_block_size, block_pointer, (current_block_size - 1) / 8 + 1, temp_sign_arr);
                block_pointer += ((current_block_size - 1) / 8 + 1);

                savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, current_block_size, temp_predict_arr, bit_count);
                block_pointer += savedbitsbytelength;
                for (j = 0; j < current_block_size; j++)
                {
                    if (temp_sign_arr[j] == 0)
                    {
                        diff = temp_predict_arr[j];
                    }
                    else
                    {
                        diff = 0 - temp_predict_arr[j];
                    }
                    current = prior + diff;
                    ori_current = (float)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(float));
                    newData_perthread++;
                }
            }
        }
        free(temp_predict_arr);
        free(temp_sign_arr);
    }

#else
    printf("Error! OpenMP not supported!\n");
#endif
}

void szp_float_decompress_openmp_threadblock_randomaccess_arg(float *newData, size_t nbEle, float absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    // *newData = (float *)malloc(sizeof(float) * nbEle);
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    unsigned int threadblocksize = 0;
    unsigned int block_size = blockSize;
#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            rcp = cmpBytes + nbThreads * sizeof(size_t);
            threadblocksize = nbEle / nbThreads;
        }
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle;
        }
        float *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;

        int prior = 0;
        int current = 0;
        int diff = 0;

        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *outputBytes_perthread = rcp + offsets[tid]; 
        unsigned char *block_pointer = outputBytes_perthread;

        float ori_prior = 0.0;
        float ori_current = 0.0;

        
        unsigned char *temp_sign_arr = (unsigned char *)malloc((block_size-1) * sizeof(unsigned char)); // 1 direct value and block_size - 1 diff. values
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc((block_size-1) * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);
            ori_prior = (float)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(float)); 
            newData_perthread ++;

            if (current_block_size > 1)
            {
                bit_count = block_pointer[0];
                block_pointer++;

                if (bit_count == 0)
                {
                    for (j = 0; j < current_block_size - 1; j++)
                    {
                        memcpy(newData_perthread, &ori_prior, sizeof(float));
                        newData_perthread++;
                    }
                }
                else
                {
                    convertByteArray2IntArray_fast_1b_args(current_block_size - 1, block_pointer, (current_block_size - 2) / 8 + 1, temp_sign_arr);
                    block_pointer += ((current_block_size - 2) / 8 + 1);

                    savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, current_block_size - 1, temp_predict_arr, bit_count);
                    block_pointer += savedbitsbytelength;
                    for (j = 0; j < current_block_size - 1; j++)
                    {
                        if (temp_sign_arr[j] == 0)
                        {
                            diff = temp_predict_arr[j];
                        }
                        else
                        {
                            diff = 0 - temp_predict_arr[j];
                        }
                        current = prior + diff;
                        ori_current = (float)current * absErrBound;
                        prior = current;
                        memcpy(newData_perthread, &ori_current, sizeof(float));
                        newData_perthread++;
                    }
                }
            }
        }
        free(temp_sign_arr);
        free(temp_predict_arr);
    }

#else
    printf("Error! OpenMP not supported!\n");
#endif
}
