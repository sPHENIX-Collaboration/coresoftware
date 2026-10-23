AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -DTREEUSERINFO_PREFIX=\"$(prefix)\"

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  TreeUserInfo.h

lib_LTLIBRARIES = \
  libtreeuserinfo.la

libtreeuserinfo_la_SOURCES = \
  TreeUserInfo.cc \
  TreeUserInfoRead.cc

libtreeuserinfo_la_LIBADD = \
  -lffamodules \
  -lffaobjects \
  -lfun4all \
  -lphool \
  $(ROOTLIBS)

bin_PROGRAMS = \
  extract-cdb-info

extract_cdb_info_SOURCES = extract-cdb-info.cc
extract_cdb_info_LDADD   = libtreeuserinfo.la

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtreeuserinfo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
