AUTOMAKE_OPTIONS = foreign subdir-objects

EXTRA_DIST = README.md PHGarfieldForwardMap.md PHGarfieldCompactMap.md

AM_CPPFLAGS = \
  -I$(srcdir) \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  PHGarfieldForwardTransport.h \
  PHGarfieldForwardMap.h \
  PHGarfieldAdaptiveMap.h \
  PHGarfieldCompactMap.h

noinst_HEADERS = \
  PHGarfieldCompactMapInternal.h \
  PHGarfieldFieldMapValidation.h \
  PHGarfieldFingerprint.h \
  PHGarfieldForwardSetup.h \
  PHGarfieldMapTools.h \
  PHGarfieldPrimitiveTree.h \
  PHGarfieldRootOutput.h \
  tests/PHGarfieldTestUtils.h

lib_LTLIBRARIES = libTpcForwardTransport.la libTpcForwardTransportTools.la

# Runtime readers, mandatory integrity checks and optional direct transport.
# It does not depend on the construction library below.
libTpcForwardTransport_la_SOURCES = \
  PHGarfieldForwardTransport.cc \
  PHGarfieldFieldMapValidation.cc \
  PHGarfieldForwardSetup.cc \
  PHGarfieldForwardMap.cc \
  PHGarfieldAdaptiveMap.cc \
  PHGarfieldCompactMap.cc
libTpcForwardTransport_la_LIBADD = -lPHGarfield

# Kept separate from runtime consumers. Installed for generation tools and
# downstream make-check fixtures that need synthetic adaptive construction.
libTpcForwardTransportTools_la_SOURCES = \
  tools/PHGarfieldAdaptiveBuilder.cc \
  tools/PHGarfieldAdaptiveGenerator.cc \
  tools/PHGarfieldForwardMapGeneration.cc \
  tools/PHGarfieldCompactMapPacking.cc \
  tools/PHGarfieldMapTools.cc
libTpcForwardTransportTools_la_LIBADD = libTpcForwardTransport.la

bin_PROGRAMS = \
  PHGarfieldForwardMapGenerator \
  PHGarfieldCompactMapMerger
PHGarfieldForwardMapGenerator_SOURCES = tools/PHGarfieldForwardMapGenerator.cc
PHGarfieldForwardMapGenerator_LDADD = libTpcForwardTransportTools.la libTpcForwardTransport.la -lffamodules
PHGarfieldCompactMapMerger_SOURCES = tools/PHGarfieldCompactMapMerger.cc
PHGarfieldCompactMapMerger_LDADD = libTpcForwardTransportTools.la libTpcForwardTransport.la

check_PROGRAMS = \
  test_adaptive_map \
  test_compact_map \
  test_forward_field_coordinates \
  test_forward_map_formula \
  test_forward_transport_config \
  test_payload_identity
TESTS = $(check_PROGRAMS)
LDADD = libTpcForwardTransportTools.la libTpcForwardTransport.la
test_adaptive_map_SOURCES = tests/test_adaptive_map.cc
test_compact_map_SOURCES = tests/test_compact_map.cc
test_forward_field_coordinates_SOURCES = tests/test_forward_field_coordinates.cc
test_forward_map_formula_SOURCES = tests/test_forward_map_formula.cc
test_forward_transport_config_SOURCES = tests/test_forward_transport_config.cc
test_payload_identity_SOURCES = tests/test_payload_identity.cc
