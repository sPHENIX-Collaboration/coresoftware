AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = `geant4-config --cflags`

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libphg4gdml.la 

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(XERCESCROOT)/include \
  -I$(OPT_SPHENIX)/include

libphg4gdml_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(XERCESCROOT)/lib \
  -L$(OPT_SPHENIX)/lib \
  `geant4-config --libs`


libphg4gdml_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all 

libphg4gdml_la_SOURCES = \
  PHG4GDMLWrite.cc \
  PHG4GDMLWriteDefine.cc \
  PHG4GDMLWriteMaterials.cc \
  PHG4GDMLWriteSolids.cc \
  PHG4GDMLWriteSetup.cc \
  PHG4GDMLWriteParamvol.cc \
  PHG4GDMLWriteStructure.cc \
  PHG4GDMLUtility.cc


##############################################
# please add new classes in alphabetical order

pkginclude_HEADERS = \
  PHG4GDMLUtility.hh \
  PHG4GDMLConfig.hh

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4gdml 

testexternals_g4gdml_SOURCES = testexternals.cc
testexternals_g4gdml_LDADD = libphg4gdml.la


testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

