#include "TpcRunAverageConditions.h"
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <tpcconditions/TpcConditions.h>
#include <stdexcept>
#include "TpcConditionsUtils.h"

TpcRunAverageConditions::TpcRunAverageConditions()
  : SubsysReco("TpcRunAverageConditions")
{
}
int TpcRunAverageConditions::InitRun(PHCompositeNode *topNode)
{
  const std::string payload = CDBInterface::instance()->getUrl("TPC_CONDITIONS");
  if (payload.empty())
  {
    throw std::runtime_error("TPC_CONDITIONS missing for run-average PHGarfield reference");
  }
  CDBTTree tree(payload);
  tree.LoadCalibrations();
  std::vector<double> p;
  std::vector<double> t;
  std::vector<double> s;
  std::vector<double> n;
  for (const auto &entry : tree.GetUInt64EntryMap())
  {
    const auto bco = entry.second.find("gbco");
    if (bco == entry.second.end() || bco->second == 0 ||
        bco->second == std::numeric_limits<std::uint64_t>::max())
    {
      continue;
    }
    double pressure = tree.GetFloatValue(entry.first, "gas_pressure");
    double temperature = tree.GetFloatValue(entry.first, "gas_temperature");
    if (std::isfinite(temperature) && temperature < 200.)
    {
      temperature += 273.15;
    }
    p.push_back(pressure);
    t.push_back(temperature);
    for (int side = 0; side < 2; ++side)
    {
      std::vector<double> values;
      for (const auto &channel : TpcConditionsDetail::r1Channels(side))
      {
        values.push_back(tree.GetFloatValue(entry.first, channel));
      }
      (side == 0 ? s : n).push_back(TpcConditionsDetail::median(values));
    }
  }
  const auto pressure = TpcConditionsDetail::summarize(p);
  const auto temperature = TpcConditionsDetail::summarize(t);
  const auto south = TpcConditionsDetail::summarize(s);
  const auto north = TpcConditionsDetail::summarize(n);
  if (p.empty() || pressure.n != p.size() || temperature.n != p.size() || south.n != p.size() ||
      north.n != p.size() || pressure.mean <= 0 || temperature.mean <= 0 || south.mean <= 0 ||
      north.mean <= 0)
  {
    throw std::runtime_error(
        "Incomplete or invalid P/T/R1-current snapshot; PHGarfield reference not initialized");
  }
  PHNodeIterator iterator(topNode);
  auto *par = dynamic_cast<PHCompositeNode *>(iterator.findFirst("PHCompositeNode", "PAR"));
  if (!par)
  {
    par = new PHCompositeNode("PAR");
    topNode->addNode(par);
  }
  auto *conditions = findNode::getClass<TpcConditions>(topNode, "TpcConditions");
  if (!conditions)
  {
    conditions = new TpcConditions();
    par->addNode(new PHDataNode<TpcConditions>(conditions, "TpcConditions", "Data"));
  }
  conditions->set_Pressure(pressure.mean);
  conditions->set_Temperature(temperature.mean);
  conditions->set_AverageLoadSR1(south.mean);
  conditions->set_LoadSR1(south.mean);
  conditions->set_AverageLoadNR1(north.mean);
  conditions->set_LoadNR1(north.mean);
  conditions->set_ConditionsAvailable(true);
  return Fun4AllReturnCodes::EVENT_OK;
}
