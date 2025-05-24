#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/szp.h include/szp_rw.h include/szp_defines.h include/szp_BytesToolkit.h include/szp_TypeManager.h include/szp_dataCompression.h include/szp_float.h include/szpd_float.h include/szp_helper.h

lib_LTLIBRARIES=libszp.la
libszp_la_CFLAGS=-I./include
if OPENMP
libszp_la_CFLAGS+=-fopenmp
endif
libszp_la_LDFLAGS = -version-info  1:0:0
libszp_la_SOURCES=src/szp.c src/szp_rw.c src/szp_BytesToolkit.c src/szp_TypeManager.c src/szp_dataCompression.c src/szp_float.c src/szpd_float.c src/szp_helper.c

libSZ_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libszp_la_CFLAGS) -O3 $(libszp_la_LDFLAGS) -o $(lib_LTLIBRARIES)
