AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libphool.la

BUILT_SOURCES = \
  testexternals.cc

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I`root-config --incdir`

libphool_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  -lEvent

libphool_la_SOURCES = \
  PHBase_dict.cc \
  PHCompositeNode.cc \
  PHFlag.cc \
  PHIOManager.cc \
  PHMessage.cc \
  PHNode.cc \
  PHNodeIOManager.cc \
  PHNodeIntegrate.cc \
  PHNodeIterator.cc \
  PHNodeReset.cc \
  PHObject.cc \
  PHOperation.cc \
  PHRandomSeed.cc \
  PHRawOManager.cc \
  PHTimer.cc \
  PHTimeServer.cc \
  PHTimeStamp.cc \
  recoConsts.cc

pkginclude_HEADERS =  \
  getClass.h \
  PHCompositeNode.h \
  PHDataNode.h \
  PHDataNodeIterator.h \
  PHFlag.h \
  PHIODataNode.h \
  PHIOManager.h \
  PHNode.h \
  PHNodeIOManager.h \
  PHNodeIntegrate.h \
  PHNodeOperation.h \
  PHNodeReset.h \
  PHNodeIterator.h \
  PHObject.h \
  phool.h \
  phooldefs.h \
  PHOperation.h \
  PHRandomSeed.h \
  PHPointerList.h \
  PHPointerListIterator.h \
  PHRawOManager.h \
  PHTimer.h \
  PHTimeServer.h \
  PHTimeStamp.h \
  PHTypedNodeIterator.h \
  recoConsts.h

noinst_HEADERS = \
  PHBase_LinkDef.h

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc

testexternals_LDADD = \
  libphool.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


PHBase_dict.cc : \
  PHFlag.h \
  phool.h \
  PHObject.h \
  PHRandomSeed.h \
  PHTimeServer.h \
  PHTimeStamp.h \
  recoConsts.h \
  PHBase_LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f *_dict* testexternals.cc *.pcm
