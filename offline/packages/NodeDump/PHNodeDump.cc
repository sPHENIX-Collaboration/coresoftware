#include "PHNodeDump.h"
#include "DumpObject.h"

#include "DumpBbcPmtInfoContainer.h"
#include "DumpBbcVertexMap.h"
#include "DumpCaloPacketContainer.h"
#include "DumpCaloTriggerInfo.h"
#include "DumpCdbUrlSave.h"
#include "DumpCentralityInfo.h"
#include "DumpEpdGeom.h"
#include "DumpEventHeader.h"
#include "DumpFlagSave.h"
#include "DumpGl1Packet.h"
#include "DumpGl1RawHit.h"
#include "DumpGlobalVertexMap.h"
#include "DumpInttDeadMap.h"
#include "DumpInttRawHitContainer.h"
#include "DumpJetContainer.h"
#include "DumpJetMap.h"
#include "DumpMbdGeom.h"
#include "DumpMbdOut.h"
#include "DumpMbdPmtContainer.h"
#include "DumpMbdVertexMap.h"
#include "DumpMicromegasRawHitContainer.h"
#include "DumpMvtxRawEvtHeader.h"
#include "DumpMvtxRawHitContainer.h"
#include "DumpPHFieldConfig.h"
#include "DumpPHG4BlockCellGeomContainer.h"
#include "DumpPHG4BlockGeomContainer.h"
#include "DumpPHG4CellContainer.h"
#include "DumpPHG4CylinderCellContainer.h"
#include "DumpPHG4CylinderCellGeomContainer.h"
#include "DumpPHG4CylinderGeomContainer.h"
#include "DumpPHG4HitContainer.h"
#include "DumpPHG4InEvent.h"
#include "DumpPHG4ParticleSvtxMap.h"
#include "DumpPHG4ScintillatorSlatContainer.h"
#include "DumpPHG4TpcCylinderGeomContainer.h"
#include "DumpPHG4TruthInfoContainer.h"
#include "DumpPHGenIntegral.h"
#include "DumpPHHepMCGenEventMap.h"
#include "DumpParticleFlowElementContainer.h"
#include "DumpPdbParameterMap.h"
#include "DumpPdbParameterMapContainer.h"
#include "DumpRawClusterContainer.h"
#include "DumpRawTowerContainer.h"
#include "DumpRawTowerGeomContainer.h"
#include "DumpRunHeader.h"
#include "DumpSvtxPHG4ParticleMap.h"
#include "DumpSvtxTrackMap.h"
#include "DumpSvtxVertexMap.h"
#include "DumpSyncObject.h"
#include "DumpTowerBackground.h"
#include "DumpTowerInfoContainer.h"
#include "DumpTpcRawHitContainer.h"
#include "DumpTpcSeedTrackMap.h"
#include "DumpTrackSeedContainer.h"
#include "DumpTrkrClusterContainer.h"
#include "DumpTrkrClusterCrossingAssoc.h"
#include "DumpTrkrClusterHitAssoc.h"
#include "DumpTrkrHitSetContainer.h"
#include "DumpTrkrHitTruthAssoc.h"
#include "DumpVariableArray.h"

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>

#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TObject.h>

#include <iostream>
#include <string>
#include <utility>

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

int PHNodeDump::AddIgnore(const std::string &name)
{
  if (ignore.find(name) != ignore.end())
  {
    std::cout << PHWHERE << " "
              << name << "already in ignore list" << std::endl;
    return -1;
  }
  ignore.insert(name);
  return 0;
}

int PHNodeDump::Select(const std::string &name)
{
  if (exclusive.find(name) != exclusive.end())
  {
    std::cout << PHWHERE << " "
              << name << "already in exclusive list" << std::endl;
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
  EventHeader *eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (eventheader)
  {
    evtsequence = eventheader->get_EvtSequence();
  }
  return 0;
}

void PHNodeDump::perform(PHNode *node)
{
  std::map<std::string, DumpObject *>::iterator iter;
  if (node->getType() == "PHIODataNode")
  {
    std::string NodeName = node->getName();
    iter = dumpthis.find(NodeName);
    if (iter == dumpthis.end())
    {
      std::cout << "Adding Dump Object for " << NodeName << std::endl;
      AddDumpObject(NodeName, node);
      iter = dumpthis.find(NodeName);  // update iterator
    }

    if (iter != dumpthis.end())
    {
      iter->second->process_event(node);
    }
    else
    {
      //           for (iter = dumpthis.begin(); iter != dumpthis.end(); iter++)
      //             {
      //               std::cout << "registered: " << iter->second->Name() << std::endl;
      //             }
      std::cout << "Something went wrong with adding Dump Object for " << NodeName
                << ", it should exist !! Trying to create it again" << std::endl;
      AddDumpObject(NodeName, node);
    }
  }
  return;
}

int PHNodeDump::CloseOutputFiles()
{
  std::map<std::string, DumpObject *>::iterator iter;
  for (iter = dumpthis.begin(); iter != dumpthis.end(); ++iter)
  {
    iter->second->CloseOutputFile();
  }
  return 0;
}

int PHNodeDump::AddDumpObject(const std::string &NodeName, PHNode *node)
{
  DumpObject *newdump;
  if (!exclusive.empty())
  {
    if (exclusive.find(NodeName) == exclusive.end())
    {
      std::cout << "Exclusive find: Ignoring " << NodeName << std::endl;
      newdump = new DumpObject(NodeName);
      newdump->NoOutput();
      return initdump(NodeName, newdump);
    }
  }
  if (ignore.find(NodeName) != ignore.end())
  {
    std::cout << "Ignoring " << NodeName << std::endl;
    newdump = new DumpObject(NodeName);
    newdump->NoOutput();
  }
  else
  {
    if (node->getType() == "PHIODataNode")
    {
      // need a static cast since only from DST these guys are of type PHIODataNode<TObject*>
      // when created they are normally  PHIODataNode<PHObject*> but can be anything else as well
      TObject *tmp = static_cast<TObject *>((static_cast<PHIODataNode<TObject> *>(node))->getData());
      if (tmp->InheritsFrom("BbcPmtInfoContainerV1"))
      {
        newdump = new DumpBbcPmtInfoContainer(NodeName);
      }
      else if (tmp->InheritsFrom("BbcVertexMap"))
      {
        newdump = new DumpBbcVertexMap(NodeName);
      }
      else if (tmp->InheritsFrom("CaloPacketContainer"))
      {
        newdump = new DumpCaloPacketContainer(NodeName);
      }
      else if (tmp->InheritsFrom("CaloTriggerInfo"))
      {
        newdump = new DumpCaloTriggerInfo(NodeName);
      }
      else if (tmp->InheritsFrom("CdbUrlSave"))
      {
        newdump = new DumpCdbUrlSave(NodeName);
      }
      else if (tmp->InheritsFrom("CentralityInfo"))
      {
        newdump = new DumpCentralityInfo(NodeName);
      }
      else if (tmp->InheritsFrom("EpdGeom"))
      {
        newdump = new DumpEpdGeom(NodeName);
      }
      else if (tmp->InheritsFrom("EventHeader"))
      {
        newdump = new DumpEventHeader(NodeName);
      }
      else if (tmp->InheritsFrom("FlagSave"))
      {
        newdump = new DumpFlagSave(NodeName);
      }
      else if (tmp->InheritsFrom("Gl1Packet"))
      {
        newdump = new DumpGl1Packet(NodeName);
      }
      else if (tmp->InheritsFrom("Gl1RawHit"))
      {
        newdump = new DumpGl1RawHit(NodeName);
      }
      else if (tmp->InheritsFrom("GlobalVertexMap"))
      {
        newdump = new DumpGlobalVertexMap(NodeName);
      }
      else if (tmp->InheritsFrom("InttDeadMap"))
      {
        newdump = new DumpInttDeadMap(NodeName);
      }
      else if (tmp->InheritsFrom("InttRawHitContainer"))
      {
        newdump = new DumpInttRawHitContainer(NodeName);
      }
      else if (tmp->InheritsFrom("JetMap"))
      {
        newdump = new DumpJetMap(NodeName);
      }
      else if (tmp->InheritsFrom("JetContainer"))
      {
        newdump = new DumpJetContainer(NodeName);
      }
      else if (tmp->InheritsFrom("MbdGeom"))
      {
        newdump = new DumpMbdGeom(NodeName);
      }
      else if (tmp->InheritsFrom("MbdOut"))
      {
        newdump = new DumpMbdOut(NodeName);
      }
      else if (tmp->InheritsFrom("MbdPmtContainer"))
      {
        newdump = new DumpMbdPmtContainer(NodeName);
      }
      else if (tmp->InheritsFrom("MbdVertexMap"))
      {
        newdump = new DumpMbdVertexMap(NodeName);
      }
      else if (tmp->InheritsFrom("MicromegasRawHitContainer"))
      {
        newdump = new DumpMicromegasRawHitContainer(NodeName);
      }
      else if (tmp->InheritsFrom("MvtxRawEvtHeader"))
      {
        newdump = new DumpMvtxRawEvtHeader(NodeName);
      }
      else if (tmp->InheritsFrom("MvtxRawHitContainer"))
      {
        newdump = new DumpMvtxRawHitContainer(NodeName);
      }
      else if (tmp->InheritsFrom("ParticleFlowElementContainer"))
      {
        newdump = new DumpParticleFlowElementContainer(NodeName);
      }
      else if (tmp->InheritsFrom("PdbParameterMap"))
      {
        newdump = new DumpPdbParameterMap(NodeName);
      }
      else if (tmp->InheritsFrom("PdbParameterMapContainer"))
      {
        newdump = new DumpPdbParameterMapContainer(NodeName);
      }
      else if (tmp->InheritsFrom("PHFieldConfig"))
      {
        newdump = new DumpPHFieldConfig(NodeName);
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
      else if (tmp->InheritsFrom("PHG4ParticleSvtxMap"))
      {
        newdump = new DumpPHG4ParticleSvtxMap(NodeName);
      }
      else if (tmp->InheritsFrom("PHG4ScintillatorSlatContainer"))
      {
        newdump = new DumpPHG4ScintillatorSlatContainer(NodeName);
      }
      else if (tmp->InheritsFrom("PHG4TpcCylinderGeomContainer"))
      {
        newdump = new DumpPHG4TpcCylinderGeomContainer(NodeName);
      }
      else if (tmp->InheritsFrom("PHG4TruthInfoContainer"))
      {
        newdump = new DumpPHG4TruthInfoContainer(NodeName);
      }
      else if (tmp->InheritsFrom("PHGenIntegral"))
      {
        newdump = new DumpPHGenIntegral(NodeName);
      }
      else if (tmp->InheritsFrom("PHHepMCGenEventMap"))
      {
        newdump = new DumpPHHepMCGenEventMap(NodeName);
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
      else if (tmp->InheritsFrom("SvtxPHG4ParticleMap"))
      {
        newdump = new DumpSvtxPHG4ParticleMap(NodeName);
      }
      else if (tmp->InheritsFrom("SvtxTrackMap"))
      {
        newdump = new DumpSvtxTrackMap(NodeName);
      }
      else if (tmp->InheritsFrom("SvtxVertexMap"))
      {
        newdump = new DumpSvtxVertexMap(NodeName);
      }
      else if (tmp->InheritsFrom("SyncObject"))
      {
        newdump = new DumpSyncObject(NodeName);
      }
      else if (tmp->InheritsFrom("TowerBackground"))
      {
        newdump = new DumpTowerBackground(NodeName);
      }
      else if (tmp->InheritsFrom("TowerInfoContainer"))
      {
        newdump = new DumpTowerInfoContainer(NodeName);
      }
      else if (tmp->InheritsFrom("TpcRawHitContainer"))
      {
        newdump = new DumpTpcRawHitContainer(NodeName);
      }
      else if (tmp->InheritsFrom("TpcSeedTrackMap"))
      {
        newdump = new DumpTpcSeedTrackMap(NodeName);
      }
      else if (tmp->InheritsFrom("TrackSeedContainer"))
      {
        newdump = new DumpTrackSeedContainer(NodeName);
      }
      else if (tmp->InheritsFrom("TrkrClusterContainer"))
      {
        newdump = new DumpTrkrClusterContainer(NodeName);
      }
      else if (tmp->InheritsFrom("TrkrClusterCrossingAssoc"))
      {
        newdump = new DumpTrkrClusterCrossingAssoc(NodeName);
      }
      else if (tmp->InheritsFrom("TrkrClusterHitAssoc"))
      {
        newdump = new DumpTrkrClusterHitAssoc(NodeName);
      }
      else if (tmp->InheritsFrom("TrkrHitSetContainer"))
      {
        newdump = new DumpTrkrHitSetContainer(NodeName);
      }
      else if (tmp->InheritsFrom("TrkrHitTruthAssoc"))
      {
        newdump = new DumpTrkrHitTruthAssoc(NodeName);
      }
      else if (tmp->InheritsFrom("VariableArray"))
      {
        newdump = new DumpVariableArray(NodeName);
      }
      else
      {
        std::cout << "Registering Dummy for " << NodeName
                  << ", Class: " << tmp->ClassName() << std::endl;
        newdump = new DumpObject(NodeName);
      }
    }
    else
    {
      std::cout << "ignoring PHDataNode: " << NodeName << std::endl;
      newdump = new DumpObject(NodeName);
    }
  }
  newdump->PrintEvtSeq(print_evtseq);
  return initdump(NodeName, newdump);
}

int PHNodeDump::initdump(const std::string &newnode, DumpObject *dmp)
{
  dmp->SetParentNodeDump(this);
  dmp->SetOutDir(outdir);
  dmp->SetPrecision(fp_precision);
  dmp->Init();
  dumpthis[newnode] = dmp;
  return 0;
}
