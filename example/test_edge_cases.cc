/**
 * @file test_edge_cases.cc
 * @brief Test edge cases and potential bugs in SIMD implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "szp_simd.h"

void test_empty_arrays() {
    printf("\n=== Testing Empty Arrays ===\n");

    // Test quantization with 0 elements
    float input[1];
    int output[1];
    szp_quantize_float_simd(input, output, 0, 1.0f);
    printf("Empty quantization: PASS\n");

    // Test max finding with 0 elements
    unsigned int max = szp_find_max_simd(NULL, 0);
    printf("Empty max finding: %u (expected 0)\n", max);
    printf("Result: %s\n", (max == 0) ? "PASS" : "FAIL");
}

void test_single_element() {
    printf("\n=== Testing Single Element ===\n");

    float input[1] = {3.14f};
    int output[1];
    szp_quantize_float_simd(input, output, 1, 100.0f);
    printf("Single element quantization: %d (expected 314)\n", output[0]);
    printf("Result: %s\n", (output[0] == 314) ? "PASS" : "FAIL");

    unsigned int values[1] = {42};
    unsigned int max = szp_find_max_simd(values, 1);
    printf("Single element max: %u (expected 42)\n", max);
    printf("Result: %s\n", (max == 42) ? "PASS" : "FAIL");
}

void test_non_multiple_of_8() {
    printf("\n=== Testing Non-Multiple of 8 Sizes ===\n");

    // Test sizes: 1, 3, 5, 7, 9, 15, 17
    int test_sizes[] = {1, 3, 5, 7, 9, 15, 17};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);

    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];
        float* input = (float*)malloc(n * sizeof(float));
        int* output = (int*)malloc(n * sizeof(int));

        // Initialize with known values
        for (int i = 0; i < n; i++) {
            input[i] = (float)(i + 1);
        }

        szp_quantize_float_simd(input, output, n, 10.0f);

        // Verify
        int errors = 0;
        for (int i = 0; i < n; i++) {
            int expected = (int)((i + 1) * 10.0f);
            if (output[i] != expected) {
                printf("  Error at size %d, index %d: expected %d, got %d\n",
                       n, i, expected, output[i]);
                errors++;
            }
        }

        if (errors == 0) {
            printf("Size %d: PASS\n", n);
        } else {
            printf("Size %d: FAIL (%d errors)\n", n, errors);
        }

        free(input);
        free(output);
    }
}

void test_integer_overflow() {
    printf("\n=== Testing Integer Overflow Cases ===\n");

    // Test near INT_MIN
    int quantized[4] = {-2147483648, 0, 1000000, -1000000}; // INT_MIN
    unsigned char signs[4];
    unsigned int magnitudes[4];
    unsigned int max_value;
    int prior = 0;

    szp_delta_encode_simd(quantized, 4, signs, magnitudes, &max_value, &prior);

    printf("Processed values including INT_MIN\n");
    for (int i = 0; i < 4; i++) {
        printf("  [%d] sign=%u, magnitude=%u\n", i, signs[i], magnitudes[i]);
    }
    printf("Max magnitude: %u\n", max_value);
    printf("Note: INT_MIN negation may cause overflow\n");
}

void test_all_zeros() {
    printf("\n=== Testing All Zeros ===\n");

    float input[16];
    for (int i = 0; i < 16; i++) input[i] = 0.0f;

    unsigned char signs[16];
    unsigned int magnitudes[16];
    unsigned int max_value;
    int prior = 0;

    szp_quantize_and_delta_simd(input, 16, 100.0f, signs, magnitudes, &max_value, &prior);

    int errors = 0;
    for (int i = 0; i < 16; i++) {
        if (signs[i] != 0 || magnitudes[i] != 0) {
            errors++;
        }
    }

    printf("All zeros test: max=%u, errors=%d\n", max_value, errors);
    printf("Result: %s\n", (max_value == 0 && errors == 0) ? "PASS" : "FAIL");
}

void test_large_values() {
    printf("\n=== Testing Large Values ===\n");

    unsigned int values[8] = {
        UINT_MAX, UINT_MAX - 1, UINT_MAX - 2, 1000000,
        500000, 250000, 100000, UINT_MAX / 2
    };

    unsigned int max = szp_find_max_simd(values, 8);
    printf("Max of large values: %u (expected %u)\n", max, UINT_MAX);
    printf("Result: %s\n", (max == UINT_MAX) ? "PASS" : "FAIL");
}

void test_negative_floats() {
    printf("\n=== Testing Negative Float Quantization ===\n");

    float input[8] = {-1.5f, -2.3f, -3.7f, 1.2f, -0.5f, 2.8f, -4.1f, 0.0f};
    int output[8];

    szp_quantize_float_simd(input, output, 8, 10.0f);

    printf("Negative float quantization:\n");
    for (int i = 0; i < 8; i++) {
        int expected = (int)(input[i] * 10.0f);
        printf("  input=%.1f, output=%d, expected=%d %s\n",
               input[i], output[i], expected,
               (output[i] == expected) ? "✓" : "✗");
    }
}

void test_alignment() {
    printf("\n=== Testing Unaligned Memory Access ===\n");

    // Allocate extra space and use unaligned pointers
    float* buffer = (float*)malloc(32 * sizeof(float) + 16);
    float* unaligned_input = (float*)((char*)buffer + 3); // 3-byte offset
    int* out_buffer = (int*)malloc(32 * sizeof(int) + 16);
    int* unaligned_output = (int*)((char*)out_buffer + 5); // 5-byte offset

    for (int i = 0; i < 16; i++) {
        unaligned_input[i] = (float)(i + 1);
    }

    szp_quantize_float_simd(unaligned_input, unaligned_output, 16, 1.0f);

    int errors = 0;
    for (int i = 0; i < 16; i++) {
        if (unaligned_output[i] != i + 1) {
            errors++;
        }
    }

    printf("Unaligned memory access: %s (%d errors)\n",
           errors == 0 ? "PASS" : "FAIL", errors);

    free(buffer);
    free(out_buffer);
}

int main() {
    printf("========================================\n");
    printf("SIMD Edge Case and Bug Testing\n");
    printf("========================================\n");
    printf("SIMD Level: %s\n", szp_simd_level());

    test_empty_arrays();
    test_single_element();
    test_non_multiple_of_8();
    test_integer_overflow();
    test_all_zeros();
    test_large_values();
    test_negative_floats();
    test_alignment();

    printf("\n========================================\n");
    printf("Edge case testing completed!\n");
    printf("========================================\n");

    return 0;
}
