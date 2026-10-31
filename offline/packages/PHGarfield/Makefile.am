AUTOMAKE_OPTIONS = foreign subdir-objects

EXTRA_DIST = \
  PHGarfieldForwardMap.md \
  tpc_frame_geometry.csv

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

noinst_HEADERS = \
  PHGarfieldPrimitiveTree.h \
  PHGarfieldFingerprint.h \
  PHGarfieldMapFingerprint.h \
  PHGarfieldFieldMapValidation.h \
  PHGarfieldForwardSetup.h \
  PHGarfieldRootOutput.h

pkginclude_HEADERS = \
  PHGarfield.h \
  PHGarfieldAdaptiveMap.h \
  PHGarfieldForwardMap.h \
  PHGarfieldRossegger.h

lib_LTLIBRARIES = \
  libPHGarfield.la

libPHGarfield_la_LIBADD = \
  -lffamodules \
  -lffarawobjects \
  -lcdbobjects \
  -lg4detectors \
  -lEvent \
  -lphool \
  -lphfield \
  -lGarfield \
  -lSubsysReco

libPHGarfield_la_SOURCES = \
  PHGarfield.cc \
  PHGarfieldAdaptiveBuilder.cc \
  PHGarfieldAdaptiveGenerator.cc \
  PHGarfieldAdaptiveMap.cc \
  PHGarfieldFieldMapValidation.cc \
  PHGarfieldForwardMap.cc \
  PHGarfieldForwardSetup.cc \
  PHGarfieldRegularMap.cc \
  PHGarfieldRossegger.cc

bin_PROGRAMS = \
  GasModel \
  MergeGasFiles \
  PHGarfieldForwardMapGenerator





PHGarfieldForwardMapGenerator_SOURCES = PHGarfieldForwardMapGenerator.cc
PHGarfieldForwardMapGenerator_LDADD = libPHGarfield.la

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

check_PROGRAMS = \
  test_adaptive_map \
  test_field_map_validation \
  test_forward_field_coordinates \
  test_forward_map_formula \
  test_forward_setup \
  test_forward_transport_config
TESTS = $(check_PROGRAMS)

test_adaptive_map_SOURCES = tests/test_adaptive_map.cc
test_adaptive_map_LDADD = libPHGarfield.la


test_forward_map_formula_SOURCES = tests/test_forward_map_formula.cc
test_forward_map_formula_LDADD = libPHGarfield.la

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

# The forward configuration checks use no external payloads.
test_forward_transport_config_SOURCES = tests/test_forward_transport_config.cc
test_forward_transport_config_LDADD = libPHGarfield.la

test_field_map_validation_SOURCES = tests/test_field_map_validation.cc
test_field_map_validation_LDADD = libPHGarfield.la

# Input policy is independent of ROOT, CDB, and PHGarfield resources.
test_forward_setup_SOURCES = tests/test_forward_setup.cc PHGarfieldForwardSetup.cc
test_forward_setup_CXXFLAGS = $(AM_CXXFLAGS)


test_forward_field_coordinates_SOURCES = tests/test_forward_field_coordinates.cc
test_forward_field_coordinates_LDADD = libPHGarfield.la
