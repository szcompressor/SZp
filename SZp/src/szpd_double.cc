/**
 *  @file szpd_double.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>, Sheng Di <sdi1@anl.gov>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "szpd_double.h"
#include <assert.h>
#include <math.h>
#include "szp_TypeManager.h"
#include "szp_CompressionToolkit.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace szp;

double *szp_double_decompress_openmp_threadblock(size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    double *newData = (double *)malloc(sizeof(double) * nbEle);
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

        double *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;

        int prior = 0;
        int current = 0;
        int diff = 0;

        unsigned int bit_count = 0;
        unsigned char *block_pointer = rcp + offsets[tid];

        double ori_prior = 0.0;
        double ori_current = 0.0;

        if (lo < hi) { // Ensure thread has data to process
            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);

            ori_prior = (double)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(double)); 
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
                ori_prior = (double)prior * absErrBound;
                
                for (j = 0; j < current_block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(double));
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
                    ori_current = (double)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(double));
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

void szp_double_decompress_single_thread_arg(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
{

    
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 0;
    
    rcp = cmpBytes + nbThreads * sizeof(size_t);
    size_t block_size = blockSize;

    size_t lo = 0;
    size_t hi = nbEle;
    double *newData_perthread = newData;
    size_t i = 0;
    size_t j = 0;

    int prior = 0;
    int current = 0;
    int diff = 0;

    unsigned int bit_count = 0;
    unsigned char *block_pointer = rcp + offsets[0];

    double ori_prior = 0.0;
    double ori_current = 0.0;

    if (lo < hi) { // Ensure there is at least one element to decompress
        memcpy(&prior, block_pointer, sizeof(int));
        block_pointer += sizeof(unsigned int);

        ori_prior = (double)prior * absErrBound;
        memcpy(newData_perthread, &ori_prior, sizeof(double)); 
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
            ori_prior = (double)prior * absErrBound;
            
            for (j = 0; j < current_block_size; j++)
            {
                memcpy(newData_perthread, &ori_prior, sizeof(double));
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
                ori_current = (double)current * absErrBound;
                prior = current;
                memcpy(newData_perthread, &ori_current, sizeof(double));
                newData_perthread++;
            }
        }
    }
    
    free(temp_predict_arr);
    free(temp_sign_arr);
}


size_t szp_double_decompress_single_thread_arg_record(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
{
    size_t total_memaccess = 0;
    
    size_t *offsets = (size_t *)cmpBytes;
    unsigned char *rcp;
    unsigned int nbThreads = 1;
    
    rcp = cmpBytes + nbThreads * sizeof(size_t);
    size_t block_size = blockSize;

    size_t lo = 0;
    size_t hi = nbEle;
    double *newData_perthread = newData;
    size_t i = 0;
    size_t j = 0;

    int prior = 0;
    int current = 0;
    int diff = 0;

    unsigned int bit_count = 0;
    unsigned char *block_pointer = rcp + offsets[0]; 
    total_memaccess += sizeof(size_t); // Reading from offsets array

    double ori_prior = 0.0;
    double ori_current = 0.0;

    if (lo < hi) { // Ensure there is at least one element to decompress
        memcpy(&prior, block_pointer, sizeof(int));
        total_memaccess += sizeof(int) * 2; // read from block_pointer and write to prior
        block_pointer += sizeof(unsigned int);

        ori_prior = (double)prior * absErrBound;
        memcpy(newData_perthread, &ori_prior, sizeof(double)); 
        total_memaccess += sizeof(double) * 2; // read from ori_prior and write to newData_perthread
        newData_perthread += 1;
    }
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int savedbitsbytelength = 0;

    // Unified loop for all remaining data blocks
    for (i = lo + 1; i < hi; i = i + block_size)
    {
        size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
        if (current_block_size == 0) continue;

        bit_count = block_pointer[0];
        total_memaccess += sizeof(unsigned char); // Reading bit_count
        block_pointer++;
        
        if (bit_count == 0)
        {
            ori_prior = (double)prior * absErrBound;
            for (j = 0; j < current_block_size; j++)
            {
                memcpy(newData_perthread, &ori_prior, sizeof(double));
                total_memaccess += sizeof(double) * 2; // read from ori_prior and write to newData_perthread
                newData_perthread++;
            }
        }
        else
        {
            size_t sign_byte_len = (current_block_size - 1) / 8 + 1;
            convertByteArray2IntArray_fast_1b_args(current_block_size, block_pointer, sign_byte_len, temp_sign_arr);
            total_memaccess += sizeof(unsigned char) * sign_byte_len;         // Reading from block_pointer
            total_memaccess += sizeof(unsigned char) * current_block_size;    // Writing to temp_sign_arr
            block_pointer += sign_byte_len;

            savedbitsbytelength = Jiajun_extract_fixed_length_bits(block_pointer, current_block_size, temp_predict_arr, bit_count);
            total_memaccess += sizeof(unsigned char) * savedbitsbytelength;   // Reading from block_pointer
            total_memaccess += sizeof(unsigned int) * current_block_size;     // Writing to temp_predict_arr
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
                total_memaccess += sizeof(unsigned char); // Reading from temp_sign_arr
                total_memaccess += sizeof(unsigned int);  // Reading from temp_predict_arr

                current = prior + diff;
                ori_current = (double)current * absErrBound;
                prior = current;
                memcpy(newData_perthread, &ori_current, sizeof(double));
                total_memaccess += sizeof(double) * 2;     // read from ori_current and write to newData_perthread
                newData_perthread++;
            }
        }
    }
    
    free(temp_predict_arr);
    free(temp_sign_arr);
    return total_memaccess;
}

void szp_double_decompress_openmp_threadblock_arg(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
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

        double *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;

        int prior = 0;
        int current = 0;
        int diff = 0;

        
        int max = 0;
        int bit_count = 0;
        unsigned char *block_pointer = rcp + offsets[tid]; 
        double ori_prior = 0.0;
        double ori_current = 0.0;
        
        if (lo < hi) { // Ensure thread has data to process
            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);

            ori_prior = (double)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(double)); 
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
                ori_prior = (double)prior * absErrBound;
                
                for (j = 0; j < current_block_size; j++)
                {
                    memcpy(newData_perthread, &ori_prior, sizeof(double));
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
                    ori_current = (double)current * absErrBound;
                    prior = current;
                    memcpy(newData_perthread, &ori_current, sizeof(double));
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

double *szp_double_decompress_openmp_threadblock_randomaccess(size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    double *newData = (double *)malloc(sizeof(double) * nbEle);
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
            hi = nbEle;
        }
        double *newData_perthread = newData + lo;
        size_t i = 0;
        size_t j = 0;   

        int prior = 0;
        int current = 0;
        int diff = 0;

        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = rcp + offsets[tid];

        double ori_prior = 0.0;
        double ori_current = 0.0;

        unsigned char *temp_sign_arr = (unsigned char *)malloc((block_size-1) * sizeof(unsigned char)); // 1 direct value and (block_size - 1) diff. values
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc((block_size-1) * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;

        for (i = lo; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);
            ori_prior = (double)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(double)); 
            newData_perthread ++;

            if (current_block_size > 1)
            {
                bit_count = block_pointer[0];
                block_pointer++;

                if (bit_count == 0)
                {
                    for (j = 0; j < current_block_size - 1; j++)
                    {
                        memcpy(newData_perthread, &ori_prior, sizeof(double));
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
                        ori_current = (double)current * absErrBound;
                        prior = current;
                        memcpy(newData_perthread, &ori_current, sizeof(double));
                        newData_perthread++;
                    }
                }
            }
        }
        free(temp_sign_arr);
        free(temp_predict_arr);
    }
    return newData;

#else
    printf("Error! OpenMP not supported!\n");
#endif
}

void szp_double_decompress_openmp_threadblock_randomaccess_arg(double *newData, size_t nbEle, double absErrBound, int blockSize, unsigned char *cmpBytes)
{
#ifdef _OPENMP
    // *newData = (double *)malloc(sizeof(double) * nbEle);
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
        if (tid == nbThreads - 1)
        {
            hi = nbEle; // last thread may have more elements
        }
        double *newData_perthread = newData + lo;
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

        double ori_prior = 0.0;
        double ori_current = 0.0;

        
        unsigned char *temp_sign_arr = (unsigned char *)malloc((block_size - 1) * sizeof(unsigned char));
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc((block_size - 1) * sizeof(unsigned int));
        unsigned int savedbitsbytelength = 0;

        for (i = lo; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i +block_size > hi) ? (hi - i) : block_size; // handle the last block size 
            if (current_block_size == 0) continue; // no data to process in this block

            memcpy(&prior, block_pointer, sizeof(int));
            block_pointer += sizeof(unsigned int);
            ori_prior = (double)prior * absErrBound;
            memcpy(newData_perthread, &ori_prior, sizeof(double));
            newData_perthread ++;

            if (current_block_size > 1)
            {
                bit_count = block_pointer[0];
                block_pointer++;
                if (bit_count == 0)
                {
                    for (j = 0; j < current_block_size - 1; j++)
                    {
                        memcpy(newData_perthread, &ori_prior, sizeof(double));
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
                        ori_current = (double)current * absErrBound;
                        prior = current;
                        memcpy(newData_perthread, &ori_current, sizeof(double));
                        newData_perthread++;
                    }
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
