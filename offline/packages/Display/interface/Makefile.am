
AUTOMAKE_OPTIONS = foreign

# List of shared libraries to produce

noinst_LTLIBRARIES = libpheve_interface.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib `root-config --evelibs`

AM_CXXFLAGS = -std=c++0x `root-config --cflags`

INCLUDES = -I$(top_srcdir)/display -I$(top_srcdir)/modules -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include -I$(OPT_PHENIX)/include 

libpheve_interface_la_LIBADD = \
	-lboost_thread \
	-lboost_system \
	-lboost_chrono \
	-lSubsysReco	

##############################################
# please add new classes in alphabetical order

include_HEADERS = \
	PHEventDisplay.h

##############################################
# please add new classes in alphabetical order

libpheve_interface_la_SOURCES = \
	PHEventDisplay.cxx \
        PHEventDisplayDict.C

# noinst_PROGRAMS = \
#   testexternals

# testexternals_SOURCES = \
#   testexternals.C

# testexternals_LDADD = \
#   libpheve.la \
#   -lboost_thread 

# testexternals.C:
# 	echo "//*** this is a generated file. Do not commit, do not edit" > $@
# 	echo "int main()" >> $@
# 	echo "{" >> $@
# 	echo "  return 0;" >> $@
# 	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES)

PHEventDisplayDict.C : PHEventDisplay.h PHEventDisplayLinkDef.h
	rootcint -f $@ -c -D__linux__ $(DEFAULT_INCLUDES) $(INCLUDES) $^
