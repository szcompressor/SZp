#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/hZCCL.h include/hZCCL_rw.h include/hZCCL_defines.h include/hZCCL_BytesToolkit.h include/hZCCL_TypeManager.h include/hZCCL_dataCompression.h include/hZCCL_utility.h include/hZCCL_float.h include/hZCCLd_float.h include/hZCCL_helper.h

lib_LTLIBRARIES=libhZCCL.la
libhZCCL_la_CFLAGS=-I./include
if OPENMP
libhZCCL_la_CFLAGS+=-fopenmp
endif
libhZCCL_la_LDFLAGS = -version-info  1:0:0
libhZCCL_la_SOURCES=src/hZCCL.c src/hZCCL_rw.c src/hZCCL_BytesToolkit.c src/hZCCL_TypeManager.c src/hZCCL_dataCompression.c src/hZCCL_utility.c src/hZCCL_float.c src/hZCCLd_float.c src/hZCCL_helper.c

libSZ_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libhZCCL_la_CFLAGS) -O3 $(libhZCCL_la_LDFLAGS) -o $(lib_LTLIBRARIES)
