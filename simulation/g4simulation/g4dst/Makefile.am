AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include

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
  -lepd_io \
  -leventplaneinfo_io \
  -lffaobjects \
  -lffarawobjects \
  -lg4detectors_io \
  -lg4eval_io \
  -lg4intt_io \
  -lg4mbd_io \
  -lg4tracking_io \
  -lglobalvertex_io \
  -lintt_io \
  -ljetbackground_io \
  -ljetbase_io \
  -lkfparticle_sphenix_io \
  -lmbd_io \
  -lmicromegas_io \
  -lmvtx_io \
  -lparticleflow_io \
  -lphfield_io \
  -lphg4hit \
  -lphgeom_io \
  -lphhepmc_io \
  -ltpc_io \
  -ltpccalib_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lzdcinfo_io

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
