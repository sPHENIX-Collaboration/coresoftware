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
  libuserinfo.la

libuserinfo_la_SOURCES = \
  TreeUserInfo.cc

libuserinfo_la_LIBADD = \
  -lffamodules \
  -lffaobjects \
  -lfun4all \
  -lphool \
  $(ROOTLIBS)

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libuserinfo.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
