#include "Timing.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>

Timing::Timing(const std::string &name)
  : SubsysReco(name)
{
}

int Timing::InitRun(PHCompositeNode * /*topNode*/)
{
  starttime = time(nullptr);
  return Fun4AllReturnCodes::EVENT_OK;
}

int Timing::process_event(PHCompositeNode * /*topNode*/)
{
  call_counter++;
  if (call_counter >= calls)
  {
    time_t difftime = time(nullptr) - starttime;
    counter++;
    std::cout << "Count " << counter << ", seconds: " << difftime << std::endl;
    starttime = time(nullptr);
    call_counter = 0;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
