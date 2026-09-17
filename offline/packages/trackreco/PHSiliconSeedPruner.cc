#include "PHSiliconSeedPruner.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

PHSiliconSeedPruner::PHSiliconSeedPruner(const std::string& name)
  : SubsysReco(name)
{
}

PHSiliconSeedPruner::~PHSiliconSeedPruner() = default;

int PHSiliconSeedPruner::Init(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedPruner::InitRun(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedPruner::process_event(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedPruner::ResetEvent(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedPruner::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
