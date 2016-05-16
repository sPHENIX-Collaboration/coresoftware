AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir` \
  -I$(GENFIT)/include

libgenfit2exp_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  -L$(GENFIT)/lib

lib_LTLIBRARIES = \
  libgenfit2exp.la

pkginclude_HEADERS = \
  fields/Field2D.h

libgenfit2exp_la_SOURCES = \
  fields/Field2D.cc

libgenfit2exp_la_LIBADD = \
  -lgenfit2

BUILT_SOURCES = \
  testexternals.C

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libgenfit2exp.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating CINT dictionaries from class headers.
#%Dict.C: %.h %LinkDef.h
#	rootcint -f $@ -c $(CINTFLAGS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

#clean-local:
#	rm -f *Dict*

testexternals_SOURCES = testexternals.C
