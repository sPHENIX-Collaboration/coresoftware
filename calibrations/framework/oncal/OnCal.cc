#include "OnCal.h"

#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/phool.h>  // for PHWHERE

#include <iostream>

OnCal::OnCal(const std::string &Name)
  : SubsysReco(Name)
  , alldone(0)
{
}

int OnCal::process_event(PHCompositeNode * /*topNode*/)
{
  std::cout << "process_event(PHCompositeNode *topNode) not implemented by daughter class: " << Name() << std::endl;
  return -1;
}

int OnCal::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "EndOfAnalysis not implemented by subsystem!" << std::endl;
  std::cout << "Use this signal for computing your calibrations and commit." << std::endl;
  std::cout << "Dont do these operations at EndOfRun since subsystems may be feeded events from different runs." << std::endl;
  std::cout << "The number of events is the real parameter here, not the runnumber." << std::endl;
  return 0;
}

void OnCal::AddComment(const std::string &adcom)
{
  if (m_Comment.size() == 0)
  {
    m_Comment = adcom;
  }
  else
  {
    m_Comment += ":";
    m_Comment += adcom;
  }
  return;
}

int OnCal::CopyTables(const int /*FromRun*/, const int /*ToRun*/, const int /*commit*/) const
{
  std::cout << PHWHERE << " CopyTables not implemented" << std::endl
            << "this calibrator cannot copy its own tables" << std::endl;
  return -1;
}
