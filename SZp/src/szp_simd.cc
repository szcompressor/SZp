#include "szp_simd.h"
#include <string.h>
#include <math.h>

// Runtime SIMD capability detection
int szp_simd_available() {
#if defined(__AVX2__)
    return 2; // AVX2 available
#elif defined(__SSE4_1__)
    return 1; // SSE4.1 available
#elif defined(__SSE2__)
    return 1; // SSE2 available
#else
    return 0; // No SIMD
#endif
}

const char* szp_simd_level() {
#if defined(__AVX2__)
    return "AVX2";
#elif defined(__SSE4_1__)
    return "SSE4.1";
#elif defined(__SSE2__)
    return "SSE2";
#else
    return "None (scalar)";
#endif
}

// Scalar fallback for quantization
static void szp_quantize_float_scalar(const float* input, int* output, size_t count, float scale) {
    for (size_t i = 0; i < count; i++) {
        output[i] = (int)(input[i] * scale);
    }
}

void szp_quantize_float_simd(const float* input, int* output, size_t count, float scale) {
#if defined(__AVX2__)
    // AVX2 implementation: process 8 floats at a time
    const size_t simd_width = 8;
    const size_t simd_count = count / simd_width;

    __m256 scale_vec = _mm256_set1_ps(scale);

    for (size_t i = 0; i < simd_count; i++) {
        // Load 8 floats
        __m256 values = _mm256_loadu_ps(&input[i * simd_width]);

        // Multiply by scale
        __m256 scaled = _mm256_mul_ps(values, scale_vec);

        // Convert to integers
        __m256i quantized = _mm256_cvtps_epi32(scaled);

        // Store 8 integers
        _mm256_storeu_si256((__m256i*)&output[i * simd_width], quantized);
    }

    // Handle remainder with scalar code
    for (size_t i = simd_count * simd_width; i < count; i++) {
        output[i] = (int)(input[i] * scale);
    }

#elif defined(__SSE2__)
    // SSE2 implementation: process 4 floats at a time
    const size_t simd_width = 4;
    const size_t simd_count = count / simd_width;

    __m128 scale_vec = _mm_set1_ps(scale);

    for (size_t i = 0; i < simd_count; i++) {
        // Load 4 floats
        __m128 values = _mm_loadu_ps(&input[i * simd_width]);

        // Multiply by scale
        __m128 scaled = _mm_mul_ps(values, scale_vec);

        // Convert to integers
        __m128i quantized = _mm_cvtps_epi32(scaled);

        // Store 4 integers
        _mm_storeu_si128((__m128i*)&output[i * simd_width], quantized);
    }

    // Handle remainder with scalar code
    for (size_t i = simd_count * simd_width; i < count; i++) {
        output[i] = (int)(input[i] * scale);
    }

#else
    // Scalar fallback
    szp_quantize_float_scalar(input, output, count, scale);
#endif
}

// Scalar fallback for max finding
static unsigned int szp_find_max_scalar(const unsigned int* values, size_t count) {
    unsigned int max_val = 0;
    for (size_t i = 0; i < count; i++) {
        if (values[i] > max_val) {
            max_val = values[i];
        }
    }
    return max_val;
}

unsigned int szp_find_max_simd(const unsigned int* values, size_t count) {
    if (count == 0) return 0;

#if defined(__AVX2__)
    // AVX2 implementation: process 8 uint32s at a time
    const size_t simd_width = 8;
    const size_t simd_count = count / simd_width;

    __m256i max_vec = _mm256_setzero_si256();

    for (size_t i = 0; i < simd_count; i++) {
        __m256i values_vec = _mm256_loadu_si256((const __m256i*)&values[i * simd_width]);
        max_vec = _mm256_max_epu32(max_vec, values_vec);
    }

    // Horizontal maximum reduction
    __m128i max_low = _mm256_castsi256_si128(max_vec);
    __m128i max_high = _mm256_extracti128_si256(max_vec, 1);
    __m128i max_128 = _mm_max_epu32(max_low, max_high);

    // Further reduce to find max of 4 elements
    __m128i shuffled = _mm_shuffle_epi32(max_128, _MM_SHUFFLE(2, 3, 0, 1));
    max_128 = _mm_max_epu32(max_128, shuffled);
    shuffled = _mm_shuffle_epi32(max_128, _MM_SHUFFLE(1, 0, 3, 2));
    max_128 = _mm_max_epu32(max_128, shuffled);

    unsigned int max_val = _mm_cvtsi128_si32(max_128);

    // Handle remainder
    for (size_t i = simd_count * simd_width; i < count; i++) {
        if (values[i] > max_val) {
            max_val = values[i];
        }
    }

    return max_val;

#elif defined(__SSE4_1__)
    // SSE4.1 implementation: process 4 uint32s at a time
    const size_t simd_width = 4;
    const size_t simd_count = count / simd_width;

    __m128i max_vec = _mm_setzero_si128();

    for (size_t i = 0; i < simd_count; i++) {
        __m128i values_vec = _mm_loadu_si128((const __m128i*)&values[i * simd_width]);
        max_vec = _mm_max_epu32(max_vec, values_vec);
    }

    // Horizontal maximum reduction
    __m128i shuffled = _mm_shuffle_epi32(max_vec, _MM_SHUFFLE(2, 3, 0, 1));
    max_vec = _mm_max_epu32(max_vec, shuffled);
    shuffled = _mm_shuffle_epi32(max_vec, _MM_SHUFFLE(1, 0, 3, 2));
    max_vec = _mm_max_epu32(max_vec, shuffled);

    unsigned int max_val = _mm_cvtsi128_si32(max_vec);

    // Handle remainder
    for (size_t i = simd_count * simd_width; i < count; i++) {
        if (values[i] > max_val) {
            max_val = values[i];
        }
    }

    return max_val;

#else
    // Scalar fallback
    return szp_find_max_scalar(values, count);
#endif
}

// Delta encoding with sign/magnitude separation
void szp_delta_encode_simd(const int* quantized, size_t count,
                           unsigned char* signs, unsigned int* magnitudes,
                           unsigned int* max_value, int* prior) {
    // This operation has data dependencies (each element depends on previous)
    // We can't fully vectorize it, but we can vectorize the sign/magnitude extraction
    // and max finding

    int prev = *prior;
    unsigned int max_val = 0;

    for (size_t i = 0; i < count; i++) {
        int current = quantized[i];
        int diff = current - prev;
        prev = current;

        if (diff == 0) {
            signs[i] = 0;
            magnitudes[i] = 0;
        } else if (diff < 0) {
            signs[i] = 1;
            magnitudes[i] = (unsigned int)(-diff);
        } else {
            signs[i] = 0;
            magnitudes[i] = (unsigned int)diff;
        }
    }

    // Use vectorized max finding
    max_val = szp_find_max_simd(magnitudes, count);

    *prior = prev;
    *max_value = max_val;
}

// Fused quantize + delta encode for better cache locality
void szp_quantize_and_delta_simd(const float* input, size_t count, float scale,
                                 unsigned char* signs, unsigned int* magnitudes,
                                 unsigned int* max_value, int* prior) {
    // NOTE: Delta encoding has a data dependency (diff = current - prev) that prevents
    // full SIMD vectorization of the quantization and differencing steps.
    // We use scalar code for quantize+delta, but SIMD for the max finding.
    // Future optimization: Consider parallel prefix-sum techniques for partial vectorization.

    int prev = *prior;

    for (size_t i = 0; i < count; i++) {
        // Quantize
        int current = (int)(input[i] * scale);

        // Delta encode
        int diff = current - prev;
        prev = current;

        // Sign/magnitude extraction
        if (diff == 0) {
            signs[i] = 0;
            magnitudes[i] = 0;
        } else if (diff < 0) {
            signs[i] = 1;
            magnitudes[i] = (unsigned int)(-diff);
        } else {
            signs[i] = 0;
            magnitudes[i] = (unsigned int)diff;
        }
    }

    // Use SIMD for max finding (this is where SIMD helps)
    unsigned int max_val = szp_find_max_simd(magnitudes, count);

    *prior = prev;
    *max_value = max_val;
}
