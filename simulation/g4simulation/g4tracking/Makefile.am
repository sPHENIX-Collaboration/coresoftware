##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

lib_LTLIBRARIES = \
   libg4tracking_io.la

libg4tracking_io_la_LIBADD = \
  -lphool \
  -lphg4hit \
  -ltrack_io

pkginclude_HEADERS = \
  TrkrTruthTrack.h \
  EmbRecoMatch.h \
  EmbRecoMatchContainer.h \
  EmbRecoMatchContainerv1.h \
  EmbRecoMatchv1.h \
  TrkrTruthTrackContainer.h \
  TrkrTruthTrackContainerv1.h \
  TrkrTruthTrackv1.h \
  TruthClusterizerBase.h

ROOTDICTS = \
  EmbRecoMatchContainer_Dict.cc \
  EmbRecoMatchContainerv1_Dict.cc \
  EmbRecoMatch_Dict.cc \
  EmbRecoMatchv1_Dict.cc \
  TrkrTruthTrackContainer_Dict.cc \
  TrkrTruthTrackContainerv1_Dict.cc \
  TrkrTruthTrack_Dict.cc \
  TrkrTruthTrackv1_Dict.cc 

pcmdir = $(libdir)
# more elegant way to create pcm files (without listing them)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

libg4tracking_io_la_SOURCES = \
  $(ROOTDICTS) \
  EmbRecoMatchContainer.cc \
  EmbRecoMatchContainerv1.cc \
  EmbRecoMatchv1.cc \
  TrkrTruthTrack.cc \
  TrkrTruthTrackContainer.cc \
  TrkrTruthTrackContainerv1.cc \
  TruthClusterizerBase.cc \
  TrkrTruthTrackv1.cc 

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_g4tracking_io

testexternals_g4tracking_io_SOURCES = testexternals.cc
testexternals_g4tracking_io_LDADD = libg4tracking_io.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
