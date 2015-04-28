AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce
lib_LTLIBRARIES = \
    libg4field.la

# there seems to be a bug in our SL5.3 gcc 4.3, when running without -fno-inline
# an stl set in PHG4FieldsPHENIX gets scrambled. There is a warning in the code when
# this is detected, if this goes away after the next compiler upgrade, this flag can be
# removed
AM_CXXFLAGS = \
  `geant4-config --cflags` -fno-inline

INCLUDES = \
    -I$(includedir) \
    -I$(OFFLINE_MAIN)/include \
    -I`root-config --incdir` \
    -I${G4_MAIN}/include

G4LIBLIST = \
    -lG4vis_management \
    -lG4intercoms \
    -lG4materials \
    -lG4particles \
    -lG4graphics_reps \
    -lG4run \
    -lG4event \
    -lG4physicslists \
    -lG4modeling \
    -lG4digits_hits \
    -lG4interfaces \
    -lG4processes \
    -lG4tracking \
    -lG4track \
    -lG4geometry \
    -lG4global

CLHEPLIST := \
    `clhep-config --libs`

libg4field_la_LDFLAGS = \
    -L$(libdir) \
    `root-config --libs` \
    `geant4-config --libs`

# Source files:

libg4field_la_SOURCES = \
    PHG4Field_Dict.C \
    PHG4Field2D.C \
    PHG4Field3D.C \
    PHG4FieldsPHENIX.cc

pkginclude_HEADERS = \
  PHG4Field2D.h \
  PHG4Field3D.h \
  PHG4FieldsPHENIX.h

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libg4field.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

# rules to make root dictionaries
PHG4Field_Dict.C: \
  PHG4Field2D.h \
  PHG4Field3D.h \
  PHG4FieldLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

