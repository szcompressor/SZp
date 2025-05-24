#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/szp.h include/szp_rw.h include/szp_defines.h include/szp_TypeManager.h include/szp_CompressionToolkit.h include/szp_float.h include/szpd_float.h include/szp_helper.h

lib_LTLIBRARIES=libszp.la
libszp_la_CXXFLAGS=-I./include
if OPENMP
libszp_la_CXXFLAGS+=-fopenmp
endif
libszp_la_LDFLAGS = -version-info  1:2:0
libszp_la_SOURCES=src/szp.cc src/szp_rw.cc src/szp_TypeManager.cc src/szp_CompressionToolkit.cc src/szp_float.cc src/szpd_float.cc src/szp_helper.cc
