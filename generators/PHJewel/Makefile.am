AUTOMAKE_OPTIONS = foreign subdir-objects

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

# JEWEL 2.6.0 / PYTHIA 6.4.25 need legacy mode under gfortran >= 10
AM_FFLAGS = -std=legacy -fallow-argument-mismatch

# ---- JEWEL core libraries, one per medium -----------------------------------
# Every core carries its own private copy of JEWEL + PYTHIA 6 and exports
# only jw_* (jewel/exports.map: the rest is local, so nothing interposes on
# PYTHIA 6 symbols in either direction).
# PHJewel dlopen()s the one named by set_medium().
noinst_LTLIBRARIES = libjewelcommon.la
libjewelcommon_la_SOURCES = \
  jewel/jewel-lib.f \
  jewel/jewel-api.f \
  jewel/pythia6425mod-lhapdf6.f \
  jewel/meix.f \
  jewel/subtraction.C \
  jewel/jewel_capi.cc
# JEWEL's own C++ is not -Werror clean; it is built as shipped
libjewelcommon_la_CXXFLAGS = -w

# an explicit version script: libtool's -export-symbols-regex is not applied
# to objects that come in through the convenience library
JEWEL_CORE_LDFLAGS = -Wl,--version-script=$(srcdir)/jewel/exports.map -L$(OFFLINE_MAIN)/lib
JEWEL_CORE_LIBADD = libjewelcommon.la -lLHAPDF -lgfortran

lib_LTLIBRARIES = libjewel_simple.la libjewel_brick.la libjewel_vac.la libPHJewel.la

libjewel_simple_la_SOURCES = jewel/medium-simple.f
libjewel_simple_la_LDFLAGS = $(JEWEL_CORE_LDFLAGS)
libjewel_simple_la_LIBADD = $(JEWEL_CORE_LIBADD)
nodist_EXTRA_libjewel_simple_la_SOURCES = dummy.cxx

libjewel_brick_la_SOURCES = jewel/medium-brick.f
libjewel_brick_la_LDFLAGS = $(JEWEL_CORE_LDFLAGS)
libjewel_brick_la_LIBADD = $(JEWEL_CORE_LIBADD)
nodist_EXTRA_libjewel_brick_la_SOURCES = dummy.cxx

libjewel_vac_la_SOURCES = jewel/medium-vac.f
libjewel_vac_la_LDFLAGS = $(JEWEL_CORE_LDFLAGS)
libjewel_vac_la_LIBADD = $(JEWEL_CORE_LIBADD)
nodist_EXTRA_libjewel_vac_la_SOURCES = dummy.cxx

# ---- the SubsysReco ---------------------------------------------------------
pkginclude_HEADERS = \
  PHJewel.h \
  PHJewelGenTrigger.h \
  PHJewelJetTrigger.h \
  PHJewelPhotonTrigger.h

libPHJewel_la_SOURCES = \
  JewelCard.cc \
  JewelCore.cc \
  JewelRun.cc \
  PHJewel.cc \
  PHJewelJetTrigger.cc \
  PHJewelPhotonTrigger.cc

libPHJewel_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs` \
  `fastjet-config --libs`

libPHJewel_la_LIBADD = \
  -lSubsysReco \
  -lphhepmc \
  -lphool \
  -lHepMC \
  -ldl

BUILT_SOURCES = \
  testexternals.cc

testexternals_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_LDADD = \
  libPHJewel.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
