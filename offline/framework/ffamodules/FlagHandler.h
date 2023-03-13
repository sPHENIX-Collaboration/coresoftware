// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_FLAGHANDLER_H
#define FFAMODULES_FLAGHANDLER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class FlagHandler : public SubsysReco
{
 public:
  FlagHandler(const std::string &name = "FlagHandler");

  ~FlagHandler() override {}

  /** Create the Flag Node if it does not exist,
      if it exists, read back flags and copy them into recoConsts
   */
  int InitRun(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

 private:
};

#endif  // FFAMODULES_FLAGHANDLER_H
