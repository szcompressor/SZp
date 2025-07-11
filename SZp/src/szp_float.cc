/**
 *  @file szp_float.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>, Sheng Di <sdi1@anl.gov>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "szp.h"
#include "szp_float.h"
#include <assert.h>
#include <math.h>
#include "szp_TypeManager.h"
#include "szp_CompressionToolkit.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace szp;

int *
szp_float_openmp_direct_predict_quantization(float *oriData, size_t *outSize, float absErrBound,
                                             size_t nbEle, int blockSize)
{
#ifdef _OPENMP

    float *op = oriData;

    size_t i = 0;


    int *quti_arr = (int *)malloc(nbEle * sizeof(int));
    int *diff_arr = (int *)malloc(nbEle * sizeof(int));
    (*outSize) = 0;


    int nbThreads = 1;
    double inver_bound = 1;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            
            inver_bound = 1 / absErrBound;
            
        }

     
#pragma omp for schedule(static)
        for (i = 0; i < nbEle; i++)
        {
            quti_arr[i] = (op[i] + absErrBound) * inver_bound;
        }

#pragma omp single
        {
            diff_arr[0] = quti_arr[0];
        }

#pragma omp for schedule(static)
        for (i = 1; i < nbEle; i++)
        {
            diff_arr[i] = quti_arr[i] - quti_arr[i - 1];
        }
    }

    free(quti_arr);

    return diff_arr;
#else
    return NULL;
#endif
}

int *
szp_float_openmp_threadblock_predict_quantization(float *oriData, size_t *outSize, float absErrBound,
                                                  size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    float *op = oriData;

  
    int *diff_arr = (int *)malloc(nbEle * sizeof(int));
    (*outSize) = 0;
    

    int nbThreads = 1;
    double inver_bound = 1;
    int threadblocksize = 1;
    int remainder = 1;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            
        }
        size_t i = 0;
    
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;
        int prior = 0;
        int current = 0;
        prior = (op[lo] + absErrBound) * inver_bound;
        diff_arr[lo] = prior;
        for (i = lo + 1; i < hi; i++)
        {
            current = (op[i] + absErrBound) * inver_bound;
            diff_arr[i] = current - prior;
            prior = current;
        }
#pragma omp single
        {
            if (remainder != 0)
            {
                size_t remainder_lo = nbEle - remainder;
                prior = (op[remainder_lo] + absErrBound) * inver_bound;
                diff_arr[remainder_lo] = prior;
                for (i = nbEle - remainder + 1; i < nbEle; i++)
                {
                    current = (op[i] + absErrBound) * inver_bound;
                    diff_arr[i] = current - prior;
                    prior = current;
                }
            }
          
        }
    }

    return diff_arr;
#else
    return NULL;
#endif
}

unsigned char *
szp_float_openmp_threadblock(float *oriData, size_t *outSize, float absErrBound,
                             size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    float *op = oriData;

    size_t maxPreservedBufferSize = sizeof(float) * nbEle + sizeof(float);
    size_t maxPreservedBufferSize_perthread = 0;

    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    unsigned char *output = (unsigned char *)malloc(maxPreservedBufferSize);
    floatToBytes(output, absErrBound); 
    unsigned char *outputBytes = output + sizeof(float); // skip the first buffer for absErrBound
   
    (*outSize) = 0;
  

    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
        }
        size_t i = 0;
        size_t j = 0;
        maxPreservedBufferSize_perthread = (sizeof(float) * nbEle + nbThreads - 1) / nbThreads;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        
        if (lo < hi) { // Ensure thread has data to process
            prior = (op[lo]) * inver_bound;
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(unsigned int);
            outSize_perthread += sizeof(unsigned int);
        }
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
        unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
       
        for (i = lo + 1; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            max = 0;
            for (j = 0; j < current_block_size; j++)
            {
                current = (op[i + j]) * inver_bound;
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                }
                else
                {
                    if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        temp_predict_arr[j] = -diff;
                    }
                    else
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = diff;
                    }
                    if (max < temp_predict_arr[j])
                        max = temp_predict_arr[j];
                }
            }

            if (max == 0) 
            {
                block_pointer[0] = 0;
                block_pointer++;
                outSize_perthread++;
            }
            else
            {
                bit_count = (int)(log2f(max)) + 1;
                block_pointer[0] = bit_count;
                
                outSize_perthread++;
                block_pointer++;
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size, block_pointer, bit_count);
                
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
#pragma omp barrier

#pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
            }
            (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
            memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
        }
#pragma omp barrier
        memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
        
        free(outputBytes_perthread);
        free(temp_sign_arr);
        free(temp_predict_arr);
#pragma omp barrier
#pragma omp single
        {
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }
    }
    (*outSize) += sizeof(float);
    return output;
#else
    printf("Error! OpenMP not supported!\n");
    return NULL;
#endif
}

/**
 * output: the first 4 bytes are used to store absErrorBound, then followed by compressed data bytes.
 * */
void szp_float_openmp_threadblock_arg(unsigned char *output, float *oriData, size_t *outSize, float absErrBound,
                                      size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    float *op = oriData;

    
    size_t maxPreservedBufferSize = sizeof(float) + sizeof(float) * nbEle; 
    size_t maxPreservedBufferSize_perthread = 0;
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    
	unsigned char* outputBytes = output + sizeof(float);
	floatToBytes(output, absErrBound);

    (*outSize) = 0;


    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
        }
        size_t i = 0;
        size_t j = 0;
        maxPreservedBufferSize_perthread = (sizeof(float) * nbEle + nbThreads - 1) / nbThreads;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        
        if (lo < hi) { // Ensure thread has data to process
            prior = (op[lo]) * inver_bound;
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(unsigned int);
            outSize_perthread += sizeof(unsigned int);
        } // if nbThreads > nbEle, threadblocksize=0, causing no data to be written
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(block_size * sizeof(unsigned char));
        unsigned int *temp_predict_arr = (unsigned int *)malloc(block_size * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo + 1; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            max = 0;
            for (j = 0; j < current_block_size; j++)
            {
                current = (op[i + j]) * inver_bound;
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                }
                else
                {
                    if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        temp_predict_arr[j] = -diff;
                    }
                    else
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = diff;
                    }
                    if (max < temp_predict_arr[j])
                        max = temp_predict_arr[j];
                }
            }

            if (max == 0) 
            {
                block_pointer[0] = 0;
                block_pointer++;
                outSize_perthread++;
            }
            else
            {
                bit_count = (int)(log2f(max)) + 1;
                block_pointer[0] = bit_count;
                
                outSize_perthread++;
                block_pointer++;
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size, block_pointer, bit_count);
                
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
#pragma omp barrier

#pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
                
            }
            (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
            memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
            
        }
#pragma omp barrier
        memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
        
        free(outputBytes_perthread);
        free(temp_sign_arr);
        free(temp_predict_arr);
#pragma omp barrier
#pragma omp single
        {
            
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }

       
    }
    
    (*outSize) += sizeof(float);

    
#else
    printf("Error! OpenMP not supported!\n");
#endif
}

void szp_float_single_thread_arg(unsigned char *output, float *oriData, size_t *outSize, float absErrBound,
                                 size_t nbEle, int blockSize)
{

    float *op = oriData;

    unsigned char* outputBytes = output + sizeof(float);
    floatToBytes(output, absErrBound);
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    (*outSize) = 0;

    double inver_bound = 1 / absErrBound;
    unsigned int block_size = blockSize;

    int nbThreads = 1;
    real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
    (*outSize) += nbThreads * sizeof(size_t); 
    outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
    offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

    size_t maxPreservedBufferSize_perthread = sizeof(float) * nbEle;
    unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
    size_t outSize_perthread = 0;
    
    int tid = 0;
    size_t lo = 0;
    size_t hi = nbEle;

    int prior = 0;
    int current = 0;
    int diff = 0;
    unsigned int max = 0;
    unsigned int bit_count = 0;
    unsigned char *block_pointer = outputBytes_perthread;
    
    if (lo < hi) {
        prior = (op[lo]) * inver_bound;
        memcpy(block_pointer, &prior, sizeof(int));
        block_pointer += sizeof(unsigned int);
        outSize_perthread += sizeof(unsigned int);
    }
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int signbytelength = 0; 
    unsigned int savedbitsbytelength = 0;
    
    for (size_t i = lo + 1; i < hi; i = i + block_size)
    {
        size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
        if (current_block_size == 0) continue;

        max = 0;
        for (size_t j = 0; j < current_block_size; j++)
        {
            current = (op[i + j]) * inver_bound;
            diff = current - prior;
            prior = current;
            if (diff == 0)
            {
                temp_sign_arr[j] = 0;
                temp_predict_arr[j] = 0;
            }
            else
            {
                if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    temp_predict_arr[j] = -diff;
                }
                else
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = diff;
                }
                if (max < temp_predict_arr[j])
                    max = temp_predict_arr[j];
            }
        }
        if (max == 0) 
        {
            block_pointer[0] = 0;
            block_pointer++;
            outSize_perthread++;
        }
        else
        {
            bit_count = (int)(log2f(max)) + 1;
            block_pointer[0] = bit_count;
            
            outSize_perthread++;
            block_pointer++;
            signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size, block_pointer); 
            block_pointer += signbytelength;
            outSize_perthread += signbytelength;
            
            savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size, block_pointer, bit_count);
            
            block_pointer += savedbitsbytelength;
            outSize_perthread += savedbitsbytelength;
        }
    }

    outSize_perthread_arr[tid] = outSize_perthread;
    offsets_perthread_arr[0] = 0;

    (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
    memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
    
    memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
    
    free(outputBytes_perthread);
    free(temp_sign_arr);
    free(temp_predict_arr);
    
    free(outSize_perthread_arr);
    free(offsets_perthread_arr);
    
    (*outSize) += sizeof(float);
}

size_t szp_float_single_thread_arg_record(unsigned char *output, float *oriData, size_t *outSize, float absErrBound,
                                       size_t nbEle, int blockSize)
{

    size_t total_memaccess = 0;
    float *op = oriData;
    
    unsigned char* outputBytes = output + sizeof(float);
    floatToBytes(output, absErrBound);
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    (*outSize) = 0;
    total_memaccess += sizeof(size_t);

    double inver_bound = 1 / absErrBound;
    unsigned int block_size = blockSize;

    int nbThreads = 1;
    real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
    (*outSize) += nbThreads * sizeof(size_t); 
   
    outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
    offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

    size_t maxPreservedBufferSize_perthread = sizeof(float) * nbEle;
    unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
    size_t outSize_perthread = 0;
    
    int tid = 0;
    size_t lo = 0;
    size_t hi = nbEle;

    int prior = 0;
    int current = 0;
    int diff = 0;
    unsigned int max = 0;
    unsigned int bit_count = 0;
    unsigned char *block_pointer = outputBytes_perthread;

    if (lo < hi) {
        prior = (op[lo]) * inver_bound;
        total_memaccess += sizeof(float);
        memcpy(block_pointer, &prior, sizeof(int));
        total_memaccess += sizeof(int) * 2; // read and write
        block_pointer += sizeof(unsigned int);
        outSize_perthread += sizeof(unsigned int);
    }
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int signbytelength = 0; 
    unsigned int savedbitsbytelength = 0;
    
    for (size_t i = lo + 1; i < hi; i = i + block_size)
    {
        size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
        if (current_block_size == 0) continue;

        max = 0;
        for (size_t j = 0; j < current_block_size; j++)
        {
            current = (op[i + j]) * inver_bound;
            total_memaccess += sizeof(float);
            diff = current - prior;
            prior = current;
            if (diff == 0)
            {
                temp_sign_arr[j] = 0;
                temp_predict_arr[j] = 0;
            }
            else
            {
                if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    temp_predict_arr[j] = -diff;
                }
                else
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = diff;
                }
                if (max < temp_predict_arr[j])
                    max = temp_predict_arr[j];
            }
            total_memaccess += sizeof(unsigned char) + sizeof(unsigned int); // for temp arrays
        }
        if (max == 0) 
        {
            block_pointer[0] = 0;
            total_memaccess += sizeof(unsigned char);
            block_pointer++;
            outSize_perthread++;
        }
        else
        {
            bit_count = (int)(log2f(max)) + 1;
            block_pointer[0] = bit_count;
            total_memaccess += sizeof(unsigned char);
            
            outSize_perthread++;
            block_pointer++;
            signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size, block_pointer); 
            total_memaccess += sizeof(unsigned char) * current_block_size; // read temp_sign_arr
            block_pointer += signbytelength;
            total_memaccess += sizeof(unsigned char) * signbytelength; // write to block_pointer
            outSize_perthread += signbytelength;
            
            savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size, block_pointer, bit_count);
            total_memaccess += sizeof(unsigned int) * current_block_size; // read temp_predict_arr
            
            block_pointer += savedbitsbytelength;
            total_memaccess += sizeof(unsigned char) * savedbitsbytelength; // write to block_pointer
            outSize_perthread += savedbitsbytelength;
        }
    }

    outSize_perthread_arr[tid] = outSize_perthread;
    total_memaccess += sizeof(size_t);

    offsets_perthread_arr[0] = 0;
    total_memaccess += sizeof(size_t);

    (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
    total_memaccess += (sizeof(size_t) * 3);
    memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
    total_memaccess += (sizeof(unsigned char) * nbThreads * sizeof(size_t)) * 2; // read and write
    
    memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
    total_memaccess += (sizeof(unsigned char) * outSize_perthread) * 2; // read and write
    
    free(outputBytes_perthread);
    free(temp_sign_arr);
    free(temp_predict_arr);
    
    free(outSize_perthread_arr);
    free(offsets_perthread_arr);
    
    (*outSize) += sizeof(float);
    return total_memaccess;
}

unsigned char *
szp_float_openmp_threadblock_randomaccess(float *oriData, size_t *outSize, float absErrBound,
                                          size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    float *op = oriData;

    size_t maxPreservedBufferSize = sizeof(float) * nbEle + sizeof(float);
    unsigned char *output = (unsigned char *)malloc(maxPreservedBufferSize);
    floatToBytes(output, absErrBound);
    
    unsigned char* outputBytes = output + sizeof(float);
    
    size_t maxPreservedBufferSize_perthread = 0;
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    
    (*outSize) = 0;
    
    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            maxPreservedBufferSize_perthread = (sizeof(float) * nbEle + nbThreads - 1) / nbThreads;
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
        }
        size_t i = 0;
        size_t j = 0;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc((block_size - 1) * sizeof(unsigned char));
        unsigned int *temp_predict_arr = (unsigned int *)malloc((block_size - 1) * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            max = 0;
            prior = (op[i]) * inver_bound;
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(unsigned int);
            outSize_perthread += sizeof(unsigned int);

            if (current_block_size > 1)
            {
                for (j = 0; j < current_block_size - 1; j++)
                {
                    current = (op[i + j + 1]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else
                    {
                        if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            temp_predict_arr[j] = -diff;
                        }
                        else
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = diff;
                        }
                        if (max < temp_predict_arr[j])
                            max = temp_predict_arr[j];
                    }
                }
            }

            if (max == 0) 
            {
                block_pointer[0] = 0;
                block_pointer++;
                outSize_perthread++;
            }
            else
            {
                bit_count = (int)(log2f(max)) + 1;
                block_pointer[0] = bit_count;
                outSize_perthread++;
                block_pointer++;
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size - 1, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size - 1, block_pointer, bit_count);
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
#pragma omp barrier

#pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
            }
            (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
            memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
        }
#pragma omp barrier
        memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
#pragma omp barrier
        
        free(outputBytes_perthread);
        free(temp_sign_arr);
        free(temp_predict_arr);
#pragma omp single
        {
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }
    }
    
    (*outSize) += sizeof(float);
    return output;
    
#else
    printf("Error! OpenMP not supported!\n");
    return NULL;
#endif
}

void
szp_float_openmp_threadblock_randomaccess_arg(unsigned char *output, float *oriData, size_t *outSize, float absErrBound,
                                          size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    float *op = oriData;

    unsigned char* outputBytes = output + sizeof(float);
    floatToBytes(output, absErrBound);
    
    size_t maxPreservedBufferSize = sizeof(float) + sizeof(float) * nbEle; 
    size_t maxPreservedBufferSize_perthread = 0;
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    
    (*outSize) = 0;
    

    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int block_size = blockSize;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            maxPreservedBufferSize_perthread = (sizeof(float) * nbEle + nbThreads - 1) / nbThreads;
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
        }
        size_t i = 0;
        size_t j = 0;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        size_t lo = tid * threadblocksize;
        size_t hi = (tid + 1) * threadblocksize;
        if (tid == nbThreads - 1) {
            hi = nbEle; // Ensure the last thread processes all remaining elements
        }

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc((block_size - 1) * sizeof(unsigned char));
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc((block_size - 1) * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        for (i = lo; i < hi; i = i + block_size)
        {
            size_t current_block_size = (i + block_size > hi) ? (hi - i) : block_size;
            if (current_block_size == 0) continue;

            max = 0;
            prior = (op[i]) * inver_bound;
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(unsigned int);
            outSize_perthread += sizeof(unsigned int);

            if (current_block_size > 1)
            {
                for (j = 0; j < current_block_size - 1; j++)
                {
                    current = (op[i + j + 1]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else
                    {
                        if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            temp_predict_arr[j] = -diff;
                        }
                        else
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = diff;
                        }
                        if (max < temp_predict_arr[j])
                            max = temp_predict_arr[j];
                    }
                }
            }

            if (max == 0) 
            {
                block_pointer[0] = 0;
                block_pointer++;
                outSize_perthread++;
            }
            else
            {
                bit_count = (int)(log2f(max)) + 1;
                block_pointer[0] = bit_count;
                outSize_perthread++;
                block_pointer++;
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, current_block_size - 1, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, current_block_size - 1, block_pointer, bit_count);
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
        }

        outSize_perthread_arr[tid] = outSize_perthread;
#pragma omp barrier

#pragma omp single
        {
            offsets_perthread_arr[0] = 0;
            for (i = 1; i < nbThreads; i++)
            {
                offsets_perthread_arr[i] = offsets_perthread_arr[i - 1] + outSize_perthread_arr[i - 1];
                
            }
            (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
            memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
            
        }
#pragma omp barrier
        memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
#pragma omp barrier
        
        free(outputBytes_perthread);
        free(temp_sign_arr);
        free(temp_predict_arr);
#pragma omp single
        {
            
            free(outSize_perthread_arr);
            free(offsets_perthread_arr);
        }
        
    }
    
    (*outSize) += sizeof(float);
    
#else
    printf("Error! OpenMP not supported!\n");
#endif
}
