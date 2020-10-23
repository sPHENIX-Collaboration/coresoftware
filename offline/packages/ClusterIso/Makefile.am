AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

pkginclude_HEADERS = \
  ClusterIso.h

libclusteriso_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

lib_LTLIBRARIES = \
  libclusteriso.la

libclusteriso_la_SOURCES = \
  ClusterIso.cc

libclusteriso_la_LIBADD = \
  -lcalo_io \
  -lg4vertex_io \
  -lSubsysReco \
  -lCLHEP

noinst_PROGRAMS = \
  testexternals

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc


testexternals_LDADD = \
  libclusteriso.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
