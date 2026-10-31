AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DPROVENANCETREE_PREFIX=\"$(prefix)\"

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  TreeUserInfo.h

lib_LTLIBRARIES = \
  libprovenancetree.la

libprovenancetree_la_SOURCES = \
  TreeUserInfo.cc \
  TreeUserInfoRead.cc

libprovenancetree_la_LIBADD = \
  -lffamodules \
  -lffaobjects \
  -lfun4all \
  -lphool \
  $(ROOTLIBS)

# plain ROOT only, so they also run without the sPHENIX libraries
bin_PROGRAMS = \
  check_provenance \
  extract-cdb-info

check_provenance_SOURCES = check_provenance.cc
check_provenance_LDADD = $(ROOTLIBS)

extract_cdb_info_SOURCES = extract-cdb-info.cc TreeUserInfoRead.cc
extract_cdb_info_LDADD = $(ROOTLIBS)
# own flags so automake compiles TreeUserInfoRead.cc separately for the program (it is also in the
# libtool library, and automake refuses to share one object between the two)
extract_cdb_info_CPPFLAGS = $(AM_CPPFLAGS)

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libprovenancetree.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
