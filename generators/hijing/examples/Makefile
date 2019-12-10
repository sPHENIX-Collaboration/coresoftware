#
# Examples Makefile.
#
#                  M. Kirsanov 25.10.2006

SHELL = /bin/sh

-include ../config.mk


# Location of directories.
TMPDIR=tmp
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=$(LIBDIR)/archive
BINDIR=bin


all:
	@echo "Usage: make examplename (without extension, example: make example_1)"


# Create an executable for one of the normal test programs

example_1 example_2 example_3 example_4 example_5 example_6 example_7 \
example_8 : ../$(LIBDIRARCH)/libhijing.a ../$(LIBDIRARCH)/libhijing_dummy.a
	@mkdir -p ../$(BINDIR)
	$(FC) $(FFLAGS) -I../$(INCDIR) $@.f -c -o ../$(TMPDIR)/$@.o
	$(CXX) $(CXXFLAGS) -I../$(INCDIR) ../$(TMPDIR)/$@.o -o ../$(BINDIR)/$@.exe \
	-L../$(LIBDIRARCH) -lhijing -lhijing_dummy \
	$(FLIBS)
	@ln -fs ../$(BINDIR)/$@.exe $@.exe


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf ../$(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core*
