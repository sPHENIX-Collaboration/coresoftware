AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  MetaDataTree.h

lib_LTLIBRARIES = \
  libmetadata.la

libmetadata_la_SOURCES = \
  MetaDataTree.cc \
  MetaDataTreeRead.cc

libmetadata_la_LIBADD = \
  -lffamodules \
  -lffaobjects \
  -lfun4all \
  -lphool \
  $(ROOTLIBS)

bin_PROGRAMS = \
  cdb-frozen-tag-check \
  check_provenance \
  extract-cdb-info

# check_provenance and extract-cdb-info: plain ROOT only, so they also run without the sPHENIX libraries

check_provenance_SOURCES = check_provenance.cc
check_provenance_LDADD = $(ROOTLIBS)

extract_cdb_info_SOURCES = extract-cdb-info.cc MetaDataTreeRead.cc
extract_cdb_info_LDADD = $(ROOTLIBS)
# own flags so automake compiles MetaDataTreeRead.cc separately for the program (it is also in the
# libtool library, and automake refuses to share one object between the two)
extract_cdb_info_CPPFLAGS = $(AM_CPPFLAGS)

# asks the CDB, so needs the sPHENIX CDB client
cdb_frozen_tag_check_SOURCES = cdb-frozen-tag-check.cc
cdb_frozen_tag_check_LDADD = -lsphenixnpc -lnopayloadclient $(ROOTLIBS)

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libmetadata.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
