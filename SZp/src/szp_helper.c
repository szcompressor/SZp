/**
 *  @file szp_helper.c
 *  @author Jiajun Huang <jiajunhuang19990916@gmail.com>
 *  @date Oct, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include "szp_helper.h"


BitArray* createBitArray(size_t size) {
    BitArray* bitArray = (BitArray*)malloc(sizeof(BitArray));
    size_t numBytes = (size + 7) / 8; // Round up to nearest byte
    bitArray->array = (unsigned char*)calloc(numBytes, sizeof(unsigned char));
    bitArray->size = size;
    return bitArray;
}

void setBit(BitArray* bitArray, size_t index, int value) {
    if (index >= bitArray->size) {
        fprintf(stderr, "Index out of bounds\n");
        return;
    }
    size_t byteIndex = index / 8;
    size_t bitPosition = index % 8;
    if (value) {
        bitArray->array[byteIndex] |= (1 << bitPosition);
    } else {
        bitArray->array[byteIndex] &= ~(1 << bitPosition);
    }
}

int getBit(BitArray* bitArray, size_t index) {
    if (index >= bitArray->size) {
        fprintf(stderr, "Index out of bounds\n");
        return -1;
    }
    size_t byteIndex = index / 8;
    size_t bitPosition = index % 8;
    return (bitArray->array[byteIndex] & (1 << bitPosition)) != 0;
}

void freeBitArray(BitArray* bitArray) {
    free(bitArray->array);
    free(bitArray);
}
