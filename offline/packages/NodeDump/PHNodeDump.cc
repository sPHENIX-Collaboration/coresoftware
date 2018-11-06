#include "PHNodeDump.h"
#include "DumpObject.h"

#include "DumpBbcVertexMap.h"
#include "DumpCaloTriggerInfo.h"
#include "DumpGlobalVertexMap.h"
#include "DumpJetMap.h"
#include "DumpPdbParameterMap.h"
#include "DumpPdbParameterMapContainer.h"
#include "DumpPHG4BlockGeomContainer.h"
#include "DumpPHG4BlockCellGeomContainer.h"
#include "DumpPHG4CellContainer.h"
#include "DumpPHG4CylinderCellContainer.h"
#include "DumpPHG4CylinderCellGeomContainer.h"
#include "DumpPHG4CylinderGeomContainer.h"
#include "DumpPHG4HitContainer.h"
#include "DumpPHG4InEvent.h"
#include "DumpPHG4ScintillatorSlatContainer.h"
#include "DumpPHG4TruthInfoContainer.h"
#include "DumpRawClusterContainer.h"
#include "DumpRawTowerContainer.h"
#include "DumpRawTowerGeomContainer.h"
#include "DumpRunHeader.h"
#include "DumpSvtxClusterMap.h"
#include "DumpSvtxHitMap.h"
#include "DumpSvtxTrackMap.h"
#include "DumpSvtxVertexMap.h"
#include "DumpSyncObject.h"
#include "DumpVariableArray.h"

#include <phool/getClass.h>

#include <ffaobjects/RunHeader.h>

#include <string>

using namespace std;



PHNodeDump::PHNodeDump()
{
  runnumber = -9999;
  evtsequence = -9999;
  fp_precision = -1;
  outdir = "./";
}

PHNodeDump::~PHNodeDump()
{
  ignore.clear();
  exclusive.clear();
  while (dumpthis.begin() != dumpthis.end())
    {
      delete dumpthis.begin()->second;
      dumpthis.erase(dumpthis.begin());
    }
  return;
}

int
PHNodeDump::AddIgnore(const string &name)
{
  if (ignore.find(name) != ignore.end())
    {
      cout << PHWHERE << " "
           << name << "already in ignore list" << endl;
      return -1;
    }
  ignore.insert(name);
  return 0;
}

int
PHNodeDump::Select(const string  &name)
{
  if (exclusive.find(name) != exclusive.end())
    {
      cout << PHWHERE << " "
           << name << "already in exclusive list" << endl;
      return -1;
    }
  exclusive.insert(name);
  return 0;
}

int PHNodeDump::GetGlobalVars(PHCompositeNode *topNode)
{
  RunHeader *runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runheader)
    {
      runnumber = runheader->get_RunNumber();
    }
  return 0;
}

void PHNodeDump::perform(PHNode* node)
{
  map <string, DumpObject *>::iterator iter;
  if (node->getType() == "PHIODataNode")
    {
      string NodeName = node->getName();
      iter = dumpthis.find(NodeName);
      if (iter == dumpthis.end())
        {
          cout << "Adding Dump Object for " << NodeName << endl;
          AddDumpObject(NodeName, node);
          iter = dumpthis.find(NodeName); // update iterator
        }

      if (iter != dumpthis.end())
        {
          iter->second->process_event(node);
        }
      else
        {
          //           for (iter = dumpthis.begin(); iter != dumpthis.end(); iter++)
          //             {
          //               cout << "registered: " << iter->second->Name() << endl;
          //             }
          cout << "Something went wrong with adding Dump Object for " << NodeName
               << ", it should exist !! Trying to create it again" << endl;
          AddDumpObject(NodeName, node);
        }

    }
  return ;
}

int PHNodeDump::CloseOutputFiles()
{
  map <string, DumpObject *>::iterator iter;
  for (iter = dumpthis.begin(); iter != dumpthis.end(); iter++)
    {
      iter->second->CloseOutputFile();
    }
  return 0;
}

int PHNodeDump::AddDumpObject(const string &NodeName, PHNode *node)
{
  DumpObject *newdump;
  string newnode = NodeName;
  if (!exclusive.empty())
    {
      if (exclusive.find(NodeName) == exclusive.end())
        {
          cout << "Exclusive find: Ignoring " << NodeName << endl;
          newdump = new DumpObject(NodeName);
	  newdump->NoOutput();
          goto initdump;
        }
    }
  if (ignore.find(NodeName) != ignore.end())
    {
      cout << "Ignoring " << NodeName << endl;
      newdump = new DumpObject(NodeName);
    }
  else
    {
      if (node->getType() == "PHIODataNode")
        {
          // need a static cast since only from DST these guys are of type PHIODataNode<TObject*>
          // when created they are normally  PHIODataNode<PHObject*> but can be anything else as well
          TObject *tmp = (TObject *)(static_cast <PHIODataNode<TObject> *>(node))->getData();
          if (tmp->InheritsFrom("BbcVertexMap"))
            {
              newdump = new DumpBbcVertexMap(NodeName);
            }
          else if (tmp->InheritsFrom("CaloTriggerInfo"))
            {
              newdump = new DumpCaloTriggerInfo(NodeName);
            }
          else if (tmp->InheritsFrom("GlobalVertexMap"))
            {
              newdump = new DumpGlobalVertexMap(NodeName);
            }
          else if (tmp->InheritsFrom("JetMap"))
            {
              newdump = new DumpJetMap(NodeName);
            }
          else if (tmp->InheritsFrom("PdbParameterMap"))
            {
              newdump = new DumpPdbParameterMap(NodeName);
            }
          else if (tmp->InheritsFrom("PdbParameterMapContainer"))
            {
              newdump = new DumpPdbParameterMapContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4BlockGeomContainer"))
            {
              newdump = new DumpPHG4BlockGeomContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4BlockCellGeomContainer"))
            {
              newdump = new DumpPHG4BlockCellGeomContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4CellContainer"))
            {
              newdump = new DumpPHG4CellContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4CylinderCellContainer"))
            {
              newdump = new DumpPHG4CylinderCellContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4CylinderGeomContainer"))
            {
              newdump = new DumpPHG4CylinderGeomContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4CylinderCellGeomContainer"))
            {
              newdump = new DumpPHG4CylinderCellGeomContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4HitContainer"))
            {
              newdump = new DumpPHG4HitContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4InEvent"))
            {
              newdump = new DumpPHG4InEvent(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4ScintillatorSlatContainer"))
            {
              newdump = new DumpPHG4ScintillatorSlatContainer(NodeName);
            }
          else if (tmp->InheritsFrom("PHG4TruthInfoContainer"))
            {
              newdump = new DumpPHG4TruthInfoContainer(NodeName);
            }
          else if (tmp->InheritsFrom("RawClusterContainer"))
            {
              newdump = new DumpRawClusterContainer(NodeName);
            }
          else if (tmp->InheritsFrom("RawTowerContainer"))
            {
              newdump = new DumpRawTowerContainer(NodeName);
            }
          else if (tmp->InheritsFrom("RawTowerGeomContainer"))
            {
              newdump = new DumpRawTowerGeomContainer(NodeName);
            }
          else if (tmp->InheritsFrom("RunHeader"))
            {
              newdump = new DumpRunHeader(NodeName);
            }
          else if (tmp->InheritsFrom("SvtxClusterMap"))
            {
              newdump = new DumpSvtxClusterMap(NodeName);
            }
          else if (tmp->InheritsFrom("SvtxHitMap"))
            {
              newdump = new DumpSvtxHitMap(NodeName);
            }
          else if (tmp->InheritsFrom("SvtxTrackMap"))
            {
              newdump = new DumpSvtxTrackMap(NodeName);
            }
          else if (tmp->InheritsFrom("SvtxVertexMap"))
            {
              newdump = new DumpSvtxVertexMap(NodeName);
            }
          else if (tmp->InheritsFrom("VariableArray"))
            {
              newdump = new DumpVariableArray(NodeName);
            }
          else
            {
              cout << "Registering Dummy for " << NodeName
                   << ", Class: " << tmp->ClassName() << endl;
              newdump = new DumpObject(NodeName);
            }
        }
      else
        {
          cout << "ignoring PHDataNode: " << NodeName << endl;
          newdump = new DumpObject(NodeName);
        }
    }

 initdump:
  newdump->SetParentNodeDump(this);
  newdump->SetOutDir(outdir);
  newdump->SetPrecision(fp_precision);
  newdump->Init();
  dumpthis[newnode] = newdump;
  return 0;
}

int
PHNodeDump::SetOutDir(const string &dirname)
{
  outdir = dirname;
  return 0;
}
