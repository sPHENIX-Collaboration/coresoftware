#include "Fun4AllTriggeredInputManager.h"

#include "SingleTriggeredInput.h"

#include <fun4all/Fun4AllInputManager.h>  // for Fun4AllInputManager
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <ffaobjects/SyncObject.h>    // for SyncObject
#include <ffaobjects/SyncObjectv1.h>  // for SyncObject

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/LL1Packet.h>
#include <ffarawobjects/LL1PacketContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <algorithm>  // for std::search
#include <cassert>
#include <climits>
#include <cstdlib>
#include <iostream>  // for operator<<, basic_ostream, endl
#include <utility>   // for pair

Fun4AllTriggeredInputManager::Fun4AllTriggeredInputManager(const std::string &name, const std::string &prdfnodename, const std::string &topnodename)
  : Fun4AllInputManager(name, prdfnodename, topnodename)
  , m_SyncObject(new SyncObjectv1())
{
  Fun4AllServer *se = Fun4AllServer::instance();
  m_topNode = se->topNode(TopNodeName());
}

Fun4AllTriggeredInputManager::~Fun4AllTriggeredInputManager()
{
  if (IsOpen())
  {
    fileclose();
  }
  delete m_SyncObject;
  for (auto *iter : m_TriggeredInputVector)
  {
    if (Verbosity() > 1)
    {
      std::cout << PHWHERE << " deleting " << iter->Name() << std::endl;
    }
    delete iter;
  }
}

int Fun4AllTriggeredInputManager::run(const int /*nevents*/)
{
  m_Gl1TriggeredInput->FillPool();
  for (auto *iter : m_TriggeredInputVector)
  {
    //    std::cout << "prdf input: " << iter->Name() << std::endl;
    iter->FillPool();
    if (iter->AllDone())
    {
      return -1;
    }
  }
  m_Gl1TriggeredInput->ReadEvent();
  for (auto *iter : m_TriggeredInputVector)
  {
    //    std::cout << "prdf input: " << iter->Name() << std::endl;
    iter->ReadEvent();
    if (iter->AllDone())
    {
      return -1;
    }
  }

  if (m_RunNumber == 0)
  {
    m_RunNumber = m_Gl1TriggeredInput->RunNumber();
    SetRunNumber(m_RunNumber);
  }
  if (m_Gl1TriggeredInput->AllDone())
  {
    return -1;
  }
  EventNumber(m_Gl1TriggeredInput->EventNumber());
  m_SyncObject->EventNumber(EventNumber());
  m_SyncObject->RunNumber(m_RunNumber);
  //    std::cout << "saving event on dst" << std::endl;
  return 0;
}

int Fun4AllTriggeredInputManager::fileclose()
{
  // for (auto iter : m_TriggerInputVector)
  // {
  //   delete iter;
  // }
  //  m_TriggerInputVector.clear();
  return 0;
}

void Fun4AllTriggeredInputManager::Print(const std::string &what) const
{
  std::cout << "Triggered Inputs" << std::endl;
  if (what == "ALL" || what == "INPUT")
  {
    for (auto *iter : m_TriggeredInputVector)
    {
      std::cout << "prdf input: " << iter->Name() << std::endl;
    }
  }
  return;
}

int Fun4AllTriggeredInputManager::ResetEvent()
{
  return 0;
}

int Fun4AllTriggeredInputManager::PushBackEvents(const int /*i*/)
{
  return 0;
}

int Fun4AllTriggeredInputManager::GetSyncObject(SyncObject **mastersync)
{
  // here we copy the sync object from the current file to the
  // location pointed to by mastersync. If mastersync is a 0 pointer
  // the syncobject is cloned. If mastersync allready exists the content
  // of syncobject is copied
  if (!(*mastersync))
  {
    if (m_SyncObject)
    {
      *mastersync = dynamic_cast<SyncObject *>(m_SyncObject->CloneMe());
      assert(*mastersync);
    }
  }
  else
  {
    *(*mastersync) = *m_SyncObject;  // copy syncobject content
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

int Fun4AllTriggeredInputManager::SyncIt(const SyncObject *mastersync)
{
  if (!mastersync)
  {
    std::cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << std::endl;
    std::cout << "Most likely your first file does not contain a SyncObject and the file" << std::endl;
    std::cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << std::endl;
    std::cout << "Change your macro and use the file opened by this input manager as first input" << std::endl;
    std::cout << "and you will be okay. Fun4All will not process the current configuration" << std::endl
              << std::endl;
    return Fun4AllReturnCodes::SYNC_FAIL;
  }
  int iret = m_SyncObject->Different(mastersync);
  if (iret)
  {
    std::cout << "big problem" << std::endl;
    exit(1);
  }
  return Fun4AllReturnCodes::SYNC_OK;
}

std::string Fun4AllTriggeredInputManager::GetString(const std::string &what) const
{
  std::cout << PHWHERE << " called with " << what << " , returning empty string" << std::endl;
  return "";
}

void Fun4AllTriggeredInputManager::registerGl1TriggeredInput(SingleTriggeredInput *prdfin)
{
  m_Gl1TriggeredInput = prdfin;
  prdfin->topNode(m_topNode);
  std::cout << "registering " << prdfin->Name() << std::endl;
}

void Fun4AllTriggeredInputManager::registerTriggeredInput(SingleTriggeredInput *prdfin)
{
  if (m_Gl1TriggeredInput == nullptr)
  {
    std::cout << "You need to register the Gl1 input manager first" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  m_TriggeredInputVector.push_back(prdfin);
  prdfin->Gl1Input(m_Gl1TriggeredInput);
  prdfin->topNode(m_topNode);
  std::cout << "registering " << prdfin->Name() << std::endl;
  //  prdfin->CreateDSTNode(m_topNode);
  //  prdfin->TriggerInputManager(this);
  return;
}
