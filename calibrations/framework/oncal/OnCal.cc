#include "OnCal.h"

#include <Event/Event.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <iostream>

using namespace std;

OnCal::OnCal(const string &Name)
  : SubsysReco(Name)
  , alldone(0)
{
}

int OnCal::process_event(PHCompositeNode * /*topNode*/)
{
  cout << "process_event(PHCompositeNode *topNode) not implemented by daughter class: " << Name() << endl;
  return -1;
}

int OnCal::End(PHCompositeNode * /*topNode*/)
{
  cout << "EndOfAnalysis not implemented by subsystem!" << endl;
  cout << "Use this signal for computing your calibrations and commit." << endl;
  cout << "Dont do these operations at EndOfRun since subsystems may be feeded events from different runs." << endl;
  cout << "The number of events is the real parameter here, not the runnumber." << endl;
  return 0;
}

void OnCal::AddComment(const string &adcom)
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
  cout << PHWHERE << " CopyTables not implemented" << endl
       << "this calibrator cannot copy its own tables" << endl;
  return -1;
}
