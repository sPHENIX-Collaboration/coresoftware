AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lHepMC \
  -lCLHEP \
  `root-config --evelibs`

bin_PROGRAMS = StarLight_To_HepMc

StarLight_To_HepMc_SOURCES = StarLight_To_HepMc.cc

clean-local:
	rm -f $(BUILT_SOURCES)

