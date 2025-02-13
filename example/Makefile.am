AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../SZp/include
AM_FCFLAGS = -I../SZp/include -I../SZp
if OPENMP
AM_CFLAGS+=-fopenmp
endif


bin_PROGRAMS=testfloat_compress_fastmode1 testfloat_decompress_fastmode1

testfloat_compress_fastmode1_SOURCES=testfloat_compress_fastmode1.c
testfloat_compress_fastmode1_LDADD=../SZp/.libs/libszp.a -lm
testfloat_decompress_fastmode1_SOURCES=testfloat_decompress_fastmode1.c
testfloat_decompress_fastmode1_LDADD=../SZp/.libs/libszp.a -lm
