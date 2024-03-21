#ifndef FUN4ALL_XINGSHIFT_C
#define FUN4ALL_XINGSHIFT_C

#include <xingshiftcal/XingShiftCal.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>

#include <phool/recoConsts.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libXingShiftCal.so)

void Fun4All_xingshift(const std::string &fname = "/sphenix/lustre01/sphnxpro/commissioning/GL1/cosmics/GL1_cosmics_gl1daq-00034390-0000.prdf", int nEvents = 2000)
//void Fun4All_xingshift(const std::string &fname = "/sphenix/lustre01/sphnxpro/commissioning/GL1/beam/GL1_beam_gl1daq-00024787-0000.prdf", int nEvents = 10000)
{

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

//recoConsts *rc = recoConsts::instance();

  XingShiftCal *xingshift = new XingShiftCal();
  se->registerSubsystem(xingshift);
 
  Fun4AllInputManager *In = new Fun4AllPrdfInputManager("in");
  In->AddFile(fname);
  se->registerInputManager(In);

  se->run(nEvents);
  se->End();
  se->PrintTimer();
  delete se;
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);




}

#endif
