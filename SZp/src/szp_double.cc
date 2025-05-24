/**
 *  @file szp_double.h
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>, Sheng Di <sdi1@anl.gov>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "szp.h"
#include "szp_double.h"
#include <assert.h>
#include <math.h>
#include "szp_TypeManager.h"
#include "szp_CompressionToolkit.h"

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace szp;

int *
szp_double_openmp_direct_predict_quantization(double *oriData, size_t *outSize, double absErrBound,
                                             size_t nbEle, int blockSize)
{
#ifdef _OPENMP

    double *op = oriData;

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
szp_double_openmp_threadblock_predict_quantization(double *oriData, size_t *outSize, double absErrBound,
                                                  size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    double *op = oriData;

  
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
szp_double_openmp_threadblock(double *oriData, size_t *outSize, double absErrBound,
                             size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    double *op = oriData;

    
    size_t maxPreservedBufferSize = sizeof(double) * nbEle; 
    size_t maxPreservedBufferSize_perthread = 0;
    unsigned char *outputBytes = (unsigned char *)malloc(maxPreservedBufferSize);
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
   
    (*outSize) = 0;
  

    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            maxPreservedBufferSize_perthread = (sizeof(double) * nbEle + nbThreads - 1) / nbThreads;
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
           remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize - 1) / block_size; 
            num_remainder_in_tb = (threadblocksize - 1) % block_size;
            
        }
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        prior = (op[lo]) * inver_bound;
        
        memcpy(block_pointer, &prior, sizeof(int));
        
        block_pointer += sizeof(unsigned int);
        outSize_perthread += sizeof(unsigned int);
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
       
        unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
       
        if (num_full_block_in_tb > 0)
        {
            for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
            {
                max = 0;
                for (j = 0; j < block_size; j++)
                {
                    current = (op[i + j]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                    
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
                
            }
        }
        
        if (num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                max = 0;
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    current = (op[i + j]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }

        
        if (tid == nbThreads - 1 && remainder != 0)
        {

            unsigned int num_full_block_in_rm = (remainder - 1) / block_size; 
            unsigned int num_remainder_in_rm = (remainder - 1) % block_size;
            prior = (op[hi]) * inver_bound;
            
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(int);
            outSize_perthread += sizeof(int);
            if (num_full_block_in_rm > 0)
            {
                
                for (i = hi + 1; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    max = 0;
                    for (j = 0; j < block_size; j++)
                    {
                        current = (op[i + j]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    for (j = 0; j < num_remainder_in_rm; j++)
                    {
                        current = (op[i + j]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
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

    
    return outputBytes;
#else
    printf("Error! OpenMP not supported!\n");
    return NULL;
#endif
}

void szp_double_openmp_threadblock_arg(unsigned char *outputBytes, double *oriData, size_t *outSize, double absErrBound,
                                      size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    double *op = oriData;

    
    size_t maxPreservedBufferSize = sizeof(double) * nbEle; 
    size_t maxPreservedBufferSize_perthread = 0;
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    


    (*outSize) = 0;


    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            maxPreservedBufferSize_perthread = (sizeof(double) * nbEle + nbThreads - 1) / nbThreads;
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize - 1) / block_size; 
            num_remainder_in_tb = (threadblocksize - 1) % block_size;
            
        }
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        prior = (op[lo]) * inver_bound;
        
        memcpy(block_pointer, &prior, sizeof(int));
        
        block_pointer += sizeof(unsigned int);
        outSize_perthread += sizeof(unsigned int);
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        if (num_full_block_in_tb > 0)
        {
            for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
            {
                max = 0;
                for (j = 0; j < block_size; j++)
                {
                    current = (op[i + j]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                    
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
                
            }
        }
        
        if (num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                max = 0;
                for (j = 0; j < num_remainder_in_tb; j++)
                {
                    current = (op[i + j]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }

        
        if (tid == nbThreads - 1 && remainder != 0)
        {

            unsigned int num_full_block_in_rm = (remainder - 1) / block_size; 
            unsigned int num_remainder_in_rm = (remainder - 1) % block_size;
            prior = (op[hi]) * inver_bound;
            
            memcpy(block_pointer, &prior, sizeof(int));
            block_pointer += sizeof(int);
            outSize_perthread += sizeof(int); 
            if (num_full_block_in_rm > 0)
            {
                
                for (i = hi + 1; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    max = 0;
                    for (j = 0; j < block_size; j++)
                    {
                        current = (op[i + j]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    for (j = 0; j < num_remainder_in_rm; j++)
                    {
                        current = (op[i + j]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
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

    
#else
    printf("Error! OpenMP not supported!\n");
#endif
}

void szp_double_single_thread_arg(unsigned char *outputBytes, double *oriData, size_t *outSize, double absErrBound,
                                 size_t nbEle, int blockSize)
{

    double *op = oriData;

    size_t maxPreservedBufferSize_perthread = 0;
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    (*outSize) = 0;

    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    //unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;

    int nbThreads = 1;
    real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
    (*outSize) += nbThreads * sizeof(size_t); 
    outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
    offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

    maxPreservedBufferSize_perthread = (sizeof(double) * nbEle + nbThreads - 1) / nbThreads;
    inver_bound = 1 / absErrBound;
    threadblocksize = nbEle / nbThreads;
    //remainder = nbEle % nbThreads; //this remainder was not used in the following code.
    num_full_block_in_tb = (threadblocksize - 1) / block_size; 
    num_remainder_in_tb = (threadblocksize - 1) % block_size;

    size_t i = 0;
    unsigned int j = 0;
    //size_t k = 0;
    unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
    size_t outSize_perthread = 0;
    
    int tid = 0;
    size_t lo = tid * threadblocksize;
    size_t hi = (tid + 1) * threadblocksize;

    int prior = 0;
    int current = 0;
    int diff = 0;
    int max = 0;
    unsigned int bit_count = 0;
    unsigned char *block_pointer = outputBytes_perthread;
    prior = (op[lo]) * inver_bound;
    
    memcpy(block_pointer, &prior, sizeof(int));
    
    block_pointer += sizeof(unsigned int);
    outSize_perthread += sizeof(unsigned int);
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int signbytelength = 0; 
    unsigned int savedbitsbytelength = 0;
    
    if (num_full_block_in_tb > 0)
    {
        for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
        {
            max = 0;
            for (j = 0; j < block_size; j++)
            {
                current = (op[i + j]) * inver_bound;
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                }
                else if (diff > 0)
                {
                    temp_sign_arr[j] = 0;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                }
                else if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    diff = 0 - diff;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
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
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
            
        }
    }
    
    if (num_remainder_in_tb > 0)
    {
        for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
        {
            max = 0;
            for (j = 0; j < num_remainder_in_tb; j++)
            {
                current = (op[i + j]) * inver_bound;
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                }
                else if (diff > 0)
                {
                    temp_sign_arr[j] = 0;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                }
                else if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    diff = 0 - diff;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
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
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
            }
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
}

size_t szp_double_single_thread_arg_record(unsigned char *outputBytes, double *oriData, size_t *outSize, double absErrBound,
                                       size_t nbEle, int blockSize)
{

    size_t total_memaccess = 0;

    double *op = oriData;

    size_t maxPreservedBufferSize_perthread = 0;
    
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;

    (*outSize) = 0;
    total_memaccess += sizeof(size_t);

    double inver_bound = 0;
    size_t threadblocksize = 0;
    //unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;

    int nbThreads = 1;
    real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
    (*outSize) += nbThreads * sizeof(size_t); 
   
    outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
    offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

    maxPreservedBufferSize_perthread = (sizeof(double) * nbEle + nbThreads - 1) / nbThreads;
    inver_bound = 1 / absErrBound;
    threadblocksize = nbEle / nbThreads;
    //remainder = nbEle % nbThreads;
    num_full_block_in_tb = (threadblocksize - 1) / block_size; 
    num_remainder_in_tb = (threadblocksize - 1) % block_size;

    size_t i = 0;
    unsigned int j = 0;
    unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
    size_t outSize_perthread = 0;
    
    int tid = 0;
    size_t lo = tid * threadblocksize;
    size_t hi = (tid + 1) * threadblocksize;

    int prior = 0;
    int current = 0;
    int diff = 0;
    int max = 0;
    unsigned int bit_count = 0;
    unsigned char *block_pointer = outputBytes_perthread;
    prior = (op[lo]) * inver_bound;
    total_memaccess += sizeof(double);
    
    memcpy(block_pointer, &prior, sizeof(int));
    total_memaccess += sizeof(int);
    total_memaccess += sizeof(int);
    
    block_pointer += sizeof(unsigned int);
    outSize_perthread += sizeof(unsigned int);
    
    unsigned char *temp_sign_arr = (unsigned char *)malloc(blockSize * sizeof(unsigned char));
    
    unsigned int *temp_predict_arr = (unsigned int *)malloc(blockSize * sizeof(unsigned int));
    unsigned int signbytelength = 0; 
    unsigned int savedbitsbytelength = 0;
    
    if (num_full_block_in_tb > 0)
    {
        for (i = lo + 1; i < hi - num_remainder_in_tb; i = i + block_size)
        {
            max = 0;
            for (j = 0; j < block_size; j++)
            {
                current = (op[i + j]) * inver_bound;
                total_memaccess += sizeof(double);
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                    total_memaccess += sizeof(unsigned int);
                    total_memaccess += sizeof(unsigned int);
                }
                else if (diff > 0)
                {
                    temp_sign_arr[j] = 0;
                    total_memaccess += sizeof(unsigned int);
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                    total_memaccess += sizeof(unsigned int);
                }
                else if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    total_memaccess += sizeof(unsigned int);
                    diff = 0 - diff;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                    total_memaccess += sizeof(unsigned int);
                }
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
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, blockSize, block_pointer); 
                total_memaccess += (sizeof(unsigned int) * blockSize);
                block_pointer += signbytelength;
                total_memaccess += (sizeof(unsigned char) * signbytelength);
                outSize_perthread += signbytelength;
                
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, blockSize, block_pointer, bit_count);
                total_memaccess += (sizeof(unsigned int) * blockSize);
                
                block_pointer += savedbitsbytelength;
                total_memaccess += (sizeof(unsigned char) * savedbitsbytelength);
                outSize_perthread += savedbitsbytelength;
            }
            
        }
    }
    
    if (num_remainder_in_tb > 0)
    {
        for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
        {
            max = 0;
            for (j = 0; j < num_remainder_in_tb; j++)
            {
                current = (op[i + j]) * inver_bound;
                total_memaccess += sizeof(double);
                diff = current - prior;
                prior = current;
                if (diff == 0)
                {
                    temp_sign_arr[j] = 0;
                    temp_predict_arr[j] = 0;
                    total_memaccess += sizeof(unsigned int);
                    total_memaccess += sizeof(unsigned int);
                }
                else if (diff > 0)
                {
                    temp_sign_arr[j] = 0;
                    total_memaccess += sizeof(unsigned int);
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                    total_memaccess += sizeof(unsigned int);
                }
                else if (diff < 0)
                {
                    temp_sign_arr[j] = 1;
                    total_memaccess += sizeof(unsigned int);
                    diff = 0 - diff;
                    if (diff > max)
                    {
                        max = diff;
                    }
                    temp_predict_arr[j] = diff;
                    total_memaccess += sizeof(unsigned int);
                }
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
                signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb, block_pointer); 
                block_pointer += signbytelength;
                outSize_perthread += signbytelength;
                total_memaccess += (sizeof(unsigned int) * num_remainder_in_tb);
                total_memaccess += (sizeof(unsigned char) * signbytelength);
                savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb, block_pointer, bit_count);
                block_pointer += savedbitsbytelength;
                outSize_perthread += savedbitsbytelength;
                total_memaccess += (sizeof(unsigned int) * num_remainder_in_tb);
                total_memaccess += (sizeof(unsigned char) * savedbitsbytelength);
            }
        }
    }

    outSize_perthread_arr[tid] = outSize_perthread;
    total_memaccess += sizeof(size_t);

    offsets_perthread_arr[0] = 0;
    total_memaccess += sizeof(size_t);

    (*outSize) += offsets_perthread_arr[nbThreads - 1] + outSize_perthread_arr[nbThreads - 1];
    total_memaccess += (sizeof(size_t) * 3);
    memcpy(outputBytes, offsets_perthread_arr, nbThreads * sizeof(size_t));
    total_memaccess += (sizeof(unsigned char) * nbThreads * sizeof(size_t));
    total_memaccess += (sizeof(unsigned char) * nbThreads * sizeof(size_t));
    

    memcpy(real_outputBytes + offsets_perthread_arr[tid], outputBytes_perthread, outSize_perthread);
    total_memaccess += (sizeof(unsigned char) * outSize_perthread);
    total_memaccess += (sizeof(unsigned char) * outSize_perthread);
    
    free(outputBytes_perthread);
    free(temp_sign_arr);
    free(temp_predict_arr);

    
    free(outSize_perthread_arr);
    free(offsets_perthread_arr);
    return total_memaccess;
}

unsigned char *
szp_double_openmp_threadblock_randomaccess(double *oriData, size_t *outSize, double absErrBound,
                                          size_t nbEle, int blockSize)
{
#ifdef _OPENMP
    
    double *op = oriData;

    
    size_t maxPreservedBufferSize = sizeof(double) * nbEle; 
    size_t maxPreservedBufferSize_perthread = 0;
    unsigned char *outputBytes = (unsigned char *)malloc(maxPreservedBufferSize);
    unsigned char *real_outputBytes; 
    size_t *outSize_perthread_arr;
    size_t *offsets_perthread_arr;
    
    (*outSize) = 0;
    

    unsigned int nbThreads = 0;
    double inver_bound = 0;
    unsigned int threadblocksize = 0;
    unsigned int remainder = 0;
    unsigned int block_size = blockSize;
    unsigned int new_block_size = block_size - 1;
    unsigned int num_full_block_in_tb = 0;
    unsigned int num_remainder_in_tb = 0;

#pragma omp parallel
    {
#pragma omp single
        {
            nbThreads = omp_get_num_threads();
            real_outputBytes = outputBytes + nbThreads * sizeof(size_t);
            (*outSize) += nbThreads * sizeof(size_t); 
            outSize_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));
            offsets_perthread_arr = (size_t *)malloc(nbThreads * sizeof(size_t));

            maxPreservedBufferSize_perthread = (sizeof(double) * nbEle + nbThreads - 1) / nbThreads;
            inver_bound = 1 / absErrBound;
            threadblocksize = nbEle / nbThreads;
            remainder = nbEle % nbThreads;
            num_full_block_in_tb = (threadblocksize) / block_size; 
            num_remainder_in_tb = (threadblocksize) % block_size;
            
        }
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        unsigned char *outputBytes_perthread = (unsigned char *)malloc(maxPreservedBufferSize_perthread);
        size_t outSize_perthread = 0;
        
        int tid = omp_get_thread_num();
        int lo = tid * threadblocksize;
        int hi = (tid + 1) * threadblocksize;

        int prior = 0;
        int current = 0;
        int diff = 0;
        unsigned int max = 0;
        unsigned int bit_count = 0;
        unsigned char *block_pointer = outputBytes_perthread;
        
        unsigned char *temp_sign_arr = (unsigned char *)malloc(new_block_size * sizeof(unsigned char));
        
        unsigned int *temp_predict_arr = (unsigned int *)malloc(new_block_size * sizeof(unsigned int));
        unsigned int signbytelength = 0; 
        unsigned int savedbitsbytelength = 0;
        
        if (num_full_block_in_tb > 0)
        {
            for (i = lo; i < hi - num_remainder_in_tb; i = i + block_size)
            {
                max = 0;
                prior = (op[i]) * inver_bound;
                memcpy(block_pointer, &prior, sizeof(int));
                block_pointer += sizeof(unsigned int);
                outSize_perthread += sizeof(unsigned int);
                for (j = 0; j < new_block_size; j++)
                {
                    current = (op[i + j + 1]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, new_block_size, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, new_block_size, block_pointer, bit_count);
                    
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }
        
        if (num_remainder_in_tb > 0)
        {
            for (i = hi - num_remainder_in_tb; i < hi; i = i + block_size)
            {
                prior = (op[i]) * inver_bound;
                memcpy(block_pointer, &prior, sizeof(int));
                block_pointer += sizeof(unsigned int);
                outSize_perthread += sizeof(unsigned int);
                max = 0;
                for (j = 0; j < num_remainder_in_tb - 1; j++)
                {
                    current = (op[i + j + 1]) * inver_bound;
                    diff = current - prior;
                    prior = current;
                    if (diff == 0)
                    {
                        temp_sign_arr[j] = 0;
                        temp_predict_arr[j] = 0;
                    }
                    else if (diff > 0)
                    {
                        temp_sign_arr[j] = 0;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
                    }
                    else if (diff < 0)
                    {
                        temp_sign_arr[j] = 1;
                        diff = 0 - diff;
                        if (diff > max)
                        {
                            max = diff;
                        }
                        temp_predict_arr[j] = diff;
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
                    signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb - 1, block_pointer); 
                    block_pointer += signbytelength;
                    outSize_perthread += signbytelength;
                    savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb - 1, block_pointer, bit_count);
                    block_pointer += savedbitsbytelength;
                    outSize_perthread += savedbitsbytelength;
                }
            }
        }

        
        if (tid == nbThreads - 1 && remainder != 0)
        {

            unsigned int num_full_block_in_rm = (remainder) / block_size; 
            unsigned int num_remainder_in_rm = (remainder) % block_size;
            
            if (num_full_block_in_rm > 0)
            {
                
                for (i = hi; i < nbEle - num_remainder_in_rm; i = i + block_size)
                {
                    prior = (op[i]) * inver_bound;
                    memcpy(block_pointer, &prior, sizeof(int));
                    block_pointer += sizeof(unsigned int);
                    outSize_perthread += sizeof(unsigned int);
                    max = 0;
                    for (j = 0; j < new_block_size; j++)
                    {
                        current = (op[i + j + 1]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, new_block_size, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, new_block_size, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
            }
            if (num_remainder_in_rm > 0)
            {
                
                for (i = nbEle - num_remainder_in_rm; i < nbEle; i = i + block_size)
                {
                    max = 0;
                    prior = (op[i]) * inver_bound;
                    memcpy(block_pointer, &prior, sizeof(int));
                    block_pointer += sizeof(unsigned int);
                    outSize_perthread += sizeof(unsigned int);
                    for (j = 0; j < num_remainder_in_rm - 1; j++)
                    {
                        current = (op[i + j + 1]) * inver_bound;
                        diff = current - prior;
                        prior = current;
                        if (diff == 0)
                        {
                            temp_sign_arr[j] = 0;
                            temp_predict_arr[j] = 0;
                        }
                        else if (diff > 0)
                        {
                            temp_sign_arr[j] = 0;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
                        }
                        else if (diff < 0)
                        {
                            temp_sign_arr[j] = 1;
                            diff = 0 - diff;
                            if (diff > max)
                            {
                                max = diff;
                            }
                            temp_predict_arr[j] = diff;
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
                        signbytelength = convertIntArray2ByteArray_fast_1b_args(temp_sign_arr, num_remainder_in_tb - 1, block_pointer); 
                        block_pointer += signbytelength;
                        outSize_perthread += signbytelength;
                        savedbitsbytelength = Jiajun_save_fixed_length_bits(temp_predict_arr, num_remainder_in_tb - 1, block_pointer, bit_count);
                        block_pointer += savedbitsbytelength;
                        outSize_perthread += savedbitsbytelength;
                    }
                }
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

    

    return outputBytes;
#else
    printf("Error! OpenMP not supported!\n");
    return NULL;
#endif
}
