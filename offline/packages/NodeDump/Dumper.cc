#include "Dumper.h"
#include "PHNodeDump.h"

#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>

#include <vector>

Dumper::Dumper(const std::string &name)
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
  std::vector<std::string> DumpNodeList;
  DumpNodeList.emplace_back("RUN");
  DumpNodeList.emplace_back("PAR");
  nodedump->PrintEvtSeq(0);
  for (const auto &iter : DumpNodeList)
  {
    if (nodeiter.cd(iter))
    {
      nodeiter.forEach(*nodedump);
      nodeiter.cd();
    }
  }
  nodedump->CloseOutputFiles();
  return 0;
}

void Dumper::SetOutDir(const std::string &outdir)
{
  nodedump->SetOutDir(outdir);
  return;
}

void Dumper::SetPrecision(const int digits)
{
  nodedump->SetPrecision(digits);
  return;
}

int Dumper::AddIgnore(const std::string &name)
{
  return nodedump->AddIgnore(name);
}

int Dumper::Select(const std::string &name)
{
  return nodedump->Select(name);
}
