#include "Fun4AllPrdfOutputManager.h"

#include <fun4all/Fun4AllOutputManager.h>  // for Fun4AllOutputManager

#include <phoolraw/PHRawOManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/recoConsts.h>

#include <iostream>

//______________________________________________________
Fun4AllPrdfOutputManager::Fun4AllPrdfOutputManager(const std::string &myname, const std::string &fname)
  : Fun4AllOutputManager(myname, fname)
  , m_PrdfNode(nullptr)
  , m_PrdfOutManager(nullptr)
{
  return;
}

//______________________________________________________
Fun4AllPrdfOutputManager::~Fun4AllPrdfOutputManager()
{
  delete m_PrdfOutManager;
  return;
}

//______________________________________________________
int Fun4AllPrdfOutputManager::InitPrdfNode(PHCompositeNode *top_node, const std::string &nodeName)
{
  PHNodeIterator topIter(top_node);
  m_PrdfNode = dynamic_cast<PHCompositeNode *>(topIter.findFirst("PHCompositeNode", nodeName));
  if (m_PrdfNode)
  {
    // the m_PrdfNode already exists (Pisa Input Mgr creates one also)
    return 0;
  }

  // check name wrt default
  if (nodeName != "SIMPRDF")
  {
    std::cout << "Fun4AllPrdfOutputManager::InitPrdfNode - WARNING: nodeName is \"" << nodeName << "\". most systems expect \"SIMPRDF\" and this is most likely not going to work" << std::endl;
  }

  // create node
  m_PrdfNode = new PHCompositeNode(nodeName);
  top_node->addNode(m_PrdfNode);
  return 0;
}

//______________________________________________________
int Fun4AllPrdfOutputManager::outfileopen(const std::string &fname)
{
  if (m_PrdfOutManager)
  {
    if (Verbosity())
    {
      std::cout << "Fun4AllPrdfOutputManager::outfileopen - closing file \"" << OutFileName() << "\"" << std::endl;
    }
    delete m_PrdfOutManager;
    m_PrdfOutManager = nullptr;
  }

  OutFileName(fname);
  if (Verbosity())
  {
    std::cout << "Fun4AllPrdfOutputManager::outfileopen - writing to file \"" << OutFileName() << "\"" << std::endl;
  }

  return 0;
}

//______________________________________________________
int Fun4AllPrdfOutputManager::Write(PHCompositeNode * /*node*/)
{
  // check m_PrdfNode
  if (!m_PrdfNode)
  {
    std::cout << "Fun4AllPrdfOutputManager::Write - prdfNode not initialized" << std::endl;
    return -1;
  }

  // check m_PrdfOutManager
  if (!m_PrdfOutManager)
  {
    InitPrdfManager();
  }
  if (!m_PrdfOutManager)
  {
    std::cout << "Fun4AllPrdfOutputManager::Write - prdf manager not initialized" << std::endl;
    return -1;
  }

  // write prdfNode to prdfManager
  bool prdf_status = m_PrdfOutManager->write(m_PrdfNode);
  return prdf_status ? 0 : -1;
}

//______________________________________________________
int Fun4AllPrdfOutputManager::InitPrdfManager()
{
  if (m_PrdfOutManager)
  {
    return -1;
  }

  // retrieve run number from recoConsts
  recoConsts *rc = recoConsts::instance();

  // retrieve run number
  int run_number = -1;
  if (rc->FlagExist("RUNNUMBER"))
  {
    run_number = rc->get_IntFlag("RUNNUMBER");
  }
  // buffer length (taken from offline/framework/simReco/PrdfReco.C)
  static const int buffer_length(8 * 1024 * 1024 / 4);

  // create output manager
  m_PrdfOutManager = new PHRawOManager(OutFileName(), run_number, buffer_length);
  return 0;
}
