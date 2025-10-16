/**
 * @file test_simd.cc
 * @brief Test program to verify SIMD vectorization functionality
 * @author Auto-generated for SZp SIMD support
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "szp_simd.h"

#define TEST_SIZE 10000
#define BLOCK_SIZE 64

void print_array(const char* name, const float* arr, int count) {
    printf("%s: [", name);
    for (int i = 0; i < count && i < 10; i++) {
        printf("%.2f%s", arr[i], (i < count-1 && i < 9) ? ", " : "");
    }
    if (count > 10) printf(", ...");
    printf("]\n");
}

void test_simd_quantization() {
    printf("\n=== Testing SIMD Quantization ===\n");

    const int n = 16;
    float input[16] = {1.5f, 2.3f, 3.7f, 4.1f, 5.9f, 6.2f, 7.8f, 8.4f,
                       9.1f, 10.6f, 11.3f, 12.7f, 13.2f, 14.8f, 15.4f, 16.9f};
    int output[16];
    float scale = 100.0f;

    printf("SIMD Level: %s\n", szp_simd_level());
    printf("Input data (first 10): ");
    print_array("", input, n);
    printf("Scale factor: %.1f\n", scale);

    szp_quantize_float_simd(input, output, n, scale);

    printf("Quantized output: [");
    for (int i = 0; i < n; i++) {
        printf("%d%s", output[i], (i < n-1) ? ", " : "");
    }
    printf("]\n");

    // Verify correctness
    int errors = 0;
    for (int i = 0; i < n; i++) {
        int expected = (int)(input[i] * scale);
        if (output[i] != expected) {
            printf("ERROR at index %d: expected %d, got %d\n", i, expected, output[i]);
            errors++;
        }
    }
    printf("Verification: %s (%d errors)\n", errors == 0 ? "PASS" : "FAIL", errors);
}

void test_simd_max_finding() {
    printf("\n=== Testing SIMD Max Finding ===\n");

    const int n = 32;
    unsigned int values[32];

    // Initialize with random values
    srand(42);
    unsigned int expected_max = 0;
    for (int i = 0; i < n; i++) {
        values[i] = rand() % 1000;
        if (values[i] > expected_max) {
            expected_max = values[i];
        }
    }

    printf("Finding max of %d random values (0-999)\n", n);
    unsigned int simd_max = szp_find_max_simd(values, n);

    printf("Expected max: %u\n", expected_max);
    printf("SIMD max:     %u\n", simd_max);
    printf("Verification: %s\n", (simd_max == expected_max) ? "PASS" : "FAIL");
}

void test_quantize_and_delta() {
    printf("\n=== Testing Fused Quantize + Delta Encode ===\n");

    const int n = 8;
    float input[8] = {1.0f, 1.1f, 1.3f, 1.2f, 1.5f, 1.4f, 1.7f, 1.9f};
    unsigned char signs[8];
    unsigned int magnitudes[8];
    unsigned int max_value = 0;
    int prior = 0;
    float scale = 100.0f;

    printf("Input: ");
    print_array("", input, n);
    printf("Scale: %.1f\n", scale);

    szp_quantize_and_delta_simd(input, n, scale, signs, magnitudes, &max_value, &prior);

    printf("Signs:      [");
    for (int i = 0; i < n; i++) {
        printf("%d%s", signs[i], (i < n-1) ? ", " : "");
    }
    printf("]\n");

    printf("Magnitudes: [");
    for (int i = 0; i < n; i++) {
        printf("%u%s", magnitudes[i], (i < n-1) ? ", " : "");
    }
    printf("]\n");

    printf("Max magnitude: %u\n", max_value);
    printf("Final prior:   %d\n", prior);
}

void benchmark_quantization() {
    printf("\n=== Benchmarking Quantization ===\n");

    const int n = TEST_SIZE;
    float* input = (float*)malloc(n * sizeof(float));
    int* output = (int*)malloc(n * sizeof(int));

    // Initialize with random data
    srand(12345);
    for (int i = 0; i < n; i++) {
        input[i] = ((float)rand() / RAND_MAX) * 100.0f;
    }

    float scale = 0.01f;  // Error bound = 100

    // Warm-up
    szp_quantize_float_simd(input, output, n, scale);

    // Benchmark
    const int iterations = 10000;
    clock_t start = clock();
    for (int iter = 0; iter < iterations; iter++) {
        szp_quantize_float_simd(input, output, n, scale);
    }
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double throughput = (double)(n * iterations) / elapsed / 1e6; // Million elements/sec

    printf("Processed %d elements x %d iterations\n", n, iterations);
    printf("Time elapsed: %.3f seconds\n", elapsed);
    printf("Throughput: %.2f million elements/second\n", throughput);

    free(input);
    free(output);
}

int main(int argc, char* argv[]) {
    printf("========================================\n");
    printf("SZp SIMD Functionality Test\n");
    printf("========================================\n");
    printf("SIMD Support Level: %s\n", szp_simd_level());
    printf("SIMD Available: %s\n", szp_simd_available() ? "Yes" : "No");

    test_simd_quantization();
    test_simd_max_finding();
    test_quantize_and_delta();
    benchmark_quantization();

    printf("\n========================================\n");
    printf("All tests completed!\n");
    printf("========================================\n");

    return 0;
}
