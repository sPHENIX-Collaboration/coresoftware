AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(OPT_SPHENIX)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(OPT_SPHENIX)/lib

lib_LTLIBRARIES = libFROG.la

# this is an ugly hack to avoid a circular dependency during the build
# In order to hook our DB handling up to fun4all (disconnecting at process_event),
# the DBInterface class has to link against libfun4all,
# but fun4all's input managers use FROG to locate files
# The solution is to build a convenience libffamodules.so which is not installed
# but contains a stripped down version of the DBInterface class
# This way the test executable linking succeeds, in production libFROG will
# pick up the real libffamodules.so from the LD_LIBRARY_PATH

noinst_LTLIBRARIES = libffamodules.la
noinst_HEADERS = \
  ffamodules/DBInterface.h

libffamodules_la_SOURCES = \
  ffamodules/DBInterface.cc

libffamodules_la_LIBADD = \
  -L$(OPT_SPHENIX)/lib \
  -lodbc++

libFROG_la_LIBADD = \
  libffamodules.la

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
