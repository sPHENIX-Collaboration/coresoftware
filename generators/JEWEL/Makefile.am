AUTOMAKE_OPTIONS = foreign

AM_FFLAGS = -fno-automatic

AM_LDFLAGS = \
  -L$(OPT_SPHENIX)/lhapdf-5.9.1/lib \
  -lLHAPDF

pkgdata_DATA = \
  params_auau.example.dat \
  medium.params.dat \
  params.example.dat

bin_PROGRAMS = \
  jewel-@PACKAGE_VERSION@-vac \
  jewel-@PACKAGE_VERSION@-simple

jewel_@PACKAGE_VERSION@_vac_SOURCES = \
  jewel-@PACKAGE_VERSION@.f \
  medium-vac.f \
  pythia6425mod.f \
  meix.f

jewel_@PACKAGE_VERSION@_simple_SOURCES = \
  jewel-@PACKAGE_VERSION@.f \
  medium-simple.f \
  pythia6425mod.f \
  meix.f

clean-local:
	rm -f $(BUILT_SOURCES)

