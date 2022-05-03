AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

xclient_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib64

xclient_LDADD = \
  -lxpload \
  -lstdc++fs



bin_PROGRAMS = \
  xclient

xclient_SOURCES = xclient.cc

clean-local:
	rm -f $(BUILT_SOURCES)
