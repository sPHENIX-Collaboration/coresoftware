AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(OPT_SPHENIX)/include \
  -isystem`root-config --incdir`

lib_LTLIBRARIES = libFROG.la

libFROG_la_LIBADD = \
  -L$(OPT_SPHENIX)/lib \
  -lodbc++

pkginclude_HEADERS = \
  FROG.h

libFROG_la_SOURCES = \
  FROG.cc

bin_SCRIPTS = \
  CreateDstList.pl \
  CreateFileList.pl

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libFROG.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
