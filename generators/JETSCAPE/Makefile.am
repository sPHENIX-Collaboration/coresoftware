AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(OPT_SPHENIX)/JETSCAPE-1.0/include \
  -I/usr/include/openmpi-x86_64

bin_PROGRAMS = \
  brickTest \
  FinalStateHadrons \
  FinalStatePartons \
  LBT_brickTest \
  MUSICTest \
  PythiaBrickTest \
  SimpleValidate \
  readerTest \
  testJetScape

#  hydroFileTest 
#  hydroJetTest 
#  hydroJetTestPGun 
#  freestream-milneTest 

brickTest_SOURCES = brickTest.cc
FinalStateHadrons_SOURCES = FinalStateHadrons.cc
FinalStatePartons_SOURCES = FinalStatePartons.cc	
LBT_brickTest_SOURCES = LBT_brickTest.cc
MUSICTest_SOURCES = MUSICTest.cc
PythiaBrickTest_SOURCES = PythiaBrickTest.cc
SimpleValidate_SOURCES = SimpleValidate.cc
freestream_milneTest_SOURCES = freestream-milneTest.cc
hydroFileTest_SOURCES = hydroFileTest.cc
hydroJetTest_SOURCES = hydroJetTest.cc
hydroJetTestPGun_SOURCES = hydroJetTestPGun.cc
readerTest_SOURCES = readerTest.cc
testJetScape_SOURCES =	testJetScape.cc


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OPT_SPHENIX)/hdf5-1.10.2/lib \
  -L$(OPT_SPHENIX)/JETSCAPE-1.0/lib \
  -L/usr/lib64/openmpi/lib \
  -lboost_system \
  -lGTL \
  -lhdf5 \
  $(OPT_SPHENIX)/HepMC-3.0.0/lib/libHepMC.so \
  -lhydroFromFile \
  -liSS_lib \
  -lJetScape \
  -lJetScapeThird \
  -lmpi \
  -lmpi_cxx \
  -lmusic_lib \
  -lpythia8

