AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include

lib_LTLIBRARIES = \
   libg4dst.la

# alphabetical order please, otherwise this list will be a mess
libg4dst_la_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lcalo_io \
  -lcalotrigger_io \
  -lcentrality_io \
  -ldecayfinder_io \
  -leventplane_io \
  -lffaobjects \
  -lg4bbc_io \
  -lg4detectors_io \
  -lg4intt_io \
  -lg4jets_io \
  -lg4vertex_io \
  -lintt_io \
  -ljetbackground_io \
  -lkfparticle_sphenix_io \
  -lmicromegas_io \
  -lmvtx_io \
  -lparticleflow_io \
  -lphfield_io \
  -lphg4hit \
  -lphgeom_io \
  -lphhepmc_io \
  -ltpc_io \
  -ltrack_reco_io \
  -ltrackbase_historic_io

libg4dst_la_SOURCES = \
   g4dst.cc

g4dst.cc:
	echo "//*** this is a generated empty file. Do not commit, do not edit" > $@

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.cc \
  g4dst.cc

testexternals_SOURCES = \
  testexternals.cc

testexternals_LDADD = \
  libg4dst.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)
