AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  PHGarfield.h

lib_LTLIBRARIES = \
  libPHGarfield.la

libPHGarfield_la_LIBADD = \
  -lffamodules \
  -lffarawobjects \
  -lcdbobjects \
  -lEvent \
  -lphool \
  -lphfield \
  -lGarfield \
  -lSubsysReco

libPHGarfield_la_SOURCES = \
  PHGarfield.cc

bin_PROGRAMS = \
  GasModel \
  MergeGasFiles

MergeGasFiles_SOURCES = MergeGasFiles.cc

MergeGasFiles_LDADD = \
  -lGarfield \
  -lphool

GasModel_SOURCES = GasModel.cc

GasModel_LDADD = \
  -lGarfield

################################################
BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libPHGarfield.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
