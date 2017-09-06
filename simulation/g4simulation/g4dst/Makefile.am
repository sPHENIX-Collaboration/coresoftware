AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
   libg4dst.la

libg4dst_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphg4hit \
  -lg4detectors_io \
  -lg4hough_io \
  -lcemc_io \
  -lg4vertex_io \
  -lg4jets_io \
  -lphgeom_io \
  -lphhepmc_io \
  -lphfield_io \
  -lcalotrigger_io \
  -lg4bbc_io \
  -ljetbackground_io

libg4dst_la_SOURCES = \
   g4dst.C

g4dst.C:
	echo "//*** this is a generated empty file. Do not commit, do not edit" > $@

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.C \
  g4dst.C

testexternals_LDADD = \
  libg4dst.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
