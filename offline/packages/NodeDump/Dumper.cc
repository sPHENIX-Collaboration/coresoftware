#include "Dumper.h"
#include "PHNodeDump.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHPointerListIterator.h>

#include <vector>

using namespace std;

Dumper::Dumper(const string &name)
  : SubsysReco(name)
{
  nodedump = new PHNodeDump();
  return;
}

Dumper::~Dumper()
{
  delete nodedump;
  return;
}

int Dumper::process_event(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  nodedump->GetGlobalVars(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (dstNode)
  {
    iter.cd("DST");
    iter.forEach(*nodedump);
    iter.cd();
  }
  return 0;
}

int Dumper::End(PHCompositeNode *topNode)
{
  PHNodeIterator nodeiter(topNode);
  vector<string> DumpNodeList;
  DumpNodeList.push_back("RUN");
  DumpNodeList.push_back("PAR");
  for (vector<string>::const_iterator iter = DumpNodeList.begin();
       iter != DumpNodeList.end(); ++iter)
  {
    if (nodeiter.cd(*iter))
    {
      nodeiter.forEach(*nodedump);
      nodeiter.cd();
    }
  }
  nodedump->CloseOutputFiles();
  return 0;
}

void Dumper::SetOutDir(const string &outdir)
{
  nodedump->SetOutDir(outdir);
  return;
}

void Dumper::SetPrecision(const int digits)
{
  nodedump->SetPrecision(digits);
  return;
}

int Dumper::AddIgnore(const string &name)
{
  return nodedump->AddIgnore(name);
}

int Dumper::Select(const string &name)
{
  return nodedump->Select(name);
}
