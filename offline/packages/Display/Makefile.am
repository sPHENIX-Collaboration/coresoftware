AUTOMAKE_OPTIONS = foreign

SUBDIRS = display modules interface

lib_LTLIBRARIES = libpheve.la

libpheve_la_SOURCES =
# Dummy C++ source to cause C++ linking (see automake manual)
nodist_EXTRA_libpheve_la_SOURCES = dummy.cxx
libpheve_la_LIBADD = \
	display/libpheve_display.la \
	modules/libpheve_modules.la \
	interface/libpheve_interface.la

clean-local:
	rm -f *Dict*

