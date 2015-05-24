#include "Dumper.h"
#include "PHNodeDump.h"

#include <phool/PHPointerListIterator.h>
#include <phool/PHNodeIterator.h>

using namespace std;

Dumper::Dumper(const string &name) : SubsysReco (name)
{
  nodedump = new PHNodeDump();
  return;
}

Dumper::~Dumper()
{
  delete nodedump;
  return;
}

int 
Dumper::process_event(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  nodedump->GetGlobalVars(topNode); 
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (dstNode)
    {
      DumpCompositeNode(dstNode);
    }
  return 0;
}

int
Dumper::DumpCompositeNode(PHCompositeNode *ThisTopNode)
{
  PHNodeIterator nodeiter(ThisTopNode);
  PHPointerListIterator<PHNode> iterat(nodeiter.ls());
  PHNode *thisNode;
  while ((thisNode = iterat()))
    {
      if (thisNode->getType() == "PHCompositeNode")
	{
	  DumpCompositeNode((PHCompositeNode*)thisNode);
	}
    }
  nodeiter.forEach(*nodedump);
  return 0;
}


int 
Dumper::End(PHCompositeNode *topNode)
{

  PHNodeIterator nodeiter(topNode);
  if (nodeiter.cd("RUN"))
    {
      nodeiter.forEach(*nodedump);
    }
  nodedump->CloseOutputFiles();
  return 0;
}

void
Dumper::SetOutDir(const string &outdir)
{
  nodedump->SetOutDir(outdir);
  return;
}

void
Dumper::SetPrecision(const int digits)
{
  nodedump->SetPrecision(digits);
  return;
}

int
Dumper::AddIgnore(const string &name)
{
  return nodedump->AddIgnore(name);
}

int
Dumper::Select(const string &name)
{
  return nodedump->Select(name);
}
