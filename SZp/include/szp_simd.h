#ifndef _SZP_SIMD_H
#define _SZP_SIMD_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Check for AVX2 support at compile time
#if defined(__AVX2__)
    #define SZP_HAS_AVX2 1
    #include <immintrin.h>
#elif defined(__SSE4_1__)
    #define SZP_HAS_SSE4_1 1
    #include <smmintrin.h>
#elif defined(__SSE2__)
    #define SZP_HAS_SSE2 1
    #include <emmintrin.h>
#endif

// Runtime SIMD capability detection
int szp_simd_available();
const char* szp_simd_level();

// Vectorized quantization: converts float array to quantized int array
// Performs: output[i] = input[i] * scale
void szp_quantize_float_simd(const float* input, int* output, size_t count, float scale);

// Vectorized delta encoding with sign/magnitude separation
// Computes differences, separates signs and absolute values, finds maximum
void szp_delta_encode_simd(const int* quantized, size_t count,
                           unsigned char* signs, unsigned int* magnitudes,
                           unsigned int* max_value, int* prior);

// Vectorized max finding for a block
unsigned int szp_find_max_simd(const unsigned int* values, size_t count);

// Combined quantize + delta encode (fused operation for better performance)
void szp_quantize_and_delta_simd(const float* input, size_t count, float scale,
                                 unsigned char* signs, unsigned int* magnitudes,
                                 unsigned int* max_value, int* prior);

#ifdef __cplusplus
}
#endif

#endif // _SZP_SIMD_H
