lib_LTLIBRARIES = \
  libSeamstress.la

pkginclude_HEADERS = \
  Seamstress.h \
  Needle.h \
  Pincushion.h

libSeamstress_la_SOURCES = \
  Seamstress.cpp

libSeamstress_la_LIBADD = \
  -lpthread

libSeamstress_la_LDFLAGS = -version-info 1:0:0

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = seamstress.pc
