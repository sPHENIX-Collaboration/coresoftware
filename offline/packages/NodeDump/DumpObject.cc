#include "DumpObject.h"

#include "PHNodeDump.h"

#include <fstream>
#include <iostream>
#include <string>

DumpObject::DumpObject(const std::string &NodeName)
  : ThisName(NodeName)
{
  return;
}

int DumpObject::Init()
{
  if (ThisName != "DUMMY")
  {
    OpenOutFile();
  }
  return 0;
}

int DumpObject::OpenOutFile()
{
  if (no_output)
  {
    return 0;
  }
  std::string fname = OutDir + "/" + ThisName + ".list";
  fout = new std::ofstream(fname.c_str());
  return 0;
}

void DumpObject::Print(const char * /*what*/) const
{
  std::cout << ThisName << " did not implement Print method" << std::endl;
  return;
}

int DumpObject::process_event(PHNode *myNode)
{
  if (fout && write_run_event)
  {
    if (fp_precision > 0)
    {
      (*fout).precision(fp_precision);
    }
    *fout << std::endl
          << "RunNumber: " << myNodeDump->RunNumber();
    if (print_evtseq)
    {
      *fout << ", Event Sequence: " << myNodeDump->EvtSequence();
    }
    *fout << std::endl;
  }
  process_Node(myNode);
  //  std::cout << ThisName << " did not implement process_event method" << std::endl;
  return 0;
}

int DumpObject::process_Node(PHNode * /*myNode*/)
{
  //  std::cout << ThisName << " did not implement process_event method" << std::endl;
  return 0;
}

int DumpObject::CloseOutputFile()
{
  delete fout;
  fout = nullptr;
  return 0;
}
