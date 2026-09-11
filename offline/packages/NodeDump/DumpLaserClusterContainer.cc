#include "DumpLaserClusterContainer.h"

#include <phool/PHIODataNode.h>

#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<LaserClusterContainer>;

DumpLaserClusterContainer::DumpLaserClusterContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
}

int DumpLaserClusterContainer::process_Node(PHNode *myNode)
{
  LaserClusterContainer *container{nullptr};
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);  // NOLINT(cppcoreguidelines-pro-type-static-cast-downcast)
  if (thisNode)
  {
    container = thisNode->getData();
  }
  if (container)
  {
    *fout << "size: " << container->size() << std::endl;
    LaserClusterContainer::ConstRange range = container->getClusters();
    for (auto iter = range.first; iter != range.second; ++iter)
    {
      const TrkrDefs::cluskey key = iter->first;
      const LaserCluster *cluster = iter->second;
      if (!cluster)
      {
        continue;
      }
      *fout << std::hex << "cluster key: 0x" << key << std::dec << std::endl;
      *fout << "getFitMode: " << cluster->getFitMode() << std::endl;
      *fout << "getX: " << cluster->getX() << std::endl;
      *fout << "getY: " << cluster->getY() << std::endl;
      *fout << "getZ: " << cluster->getZ() << std::endl;
      *fout << "getLayer: " << cluster->getLayer() << std::endl;
      *fout << "getIPhi: " << cluster->getIPhi() << std::endl;
      *fout << "getIT: " << cluster->getIT() << std::endl;
      *fout << "getLayerInt: " << cluster->getLayerInt() << std::endl;
      *fout << "getIPhiInt: " << cluster->getIPhiInt() << std::endl;
      *fout << "getITInt: " << cluster->getITInt() << std::endl;
      *fout << "getAdc: " << cluster->getAdc() << std::endl;
      *fout << "getNhits: " << cluster->getNhits() << std::endl;
      *fout << "getNLayers: " << cluster->getNLayers() << std::endl;
      *fout << "getNIPhi: " << cluster->getNIPhi() << std::endl;
      *fout << "getNIT: " << cluster->getNIT() << std::endl;
      *fout << "getSDLayer: " << cluster->getSDLayer() << std::endl;
      *fout << "getSDIPhi: " << cluster->getSDIPhi() << std::endl;
      *fout << "getSDIT: " << cluster->getSDIT() << std::endl;
      *fout << "getSDWeightedLayer: " << cluster->getSDWeightedLayer() << std::endl;
      *fout << "getSDWeightedIPhi: " << cluster->getSDWeightedIPhi() << std::endl;
      *fout << "getSDWeightedIT: " << cluster->getSDWeightedIT() << std::endl;

      for (unsigned int i = 0; i < cluster->getNhits(); ++i)
      {
        const LaserClusterHitInfo &hit = cluster->getHit(i);
        *fout << std::hex << "getHit(" << i << ").hitsetkey: 0x" << hit.hitsetkey
              << ", hitkey: 0x" << hit.hitkey << std::dec << ", adc: " << hit.adc << std::endl;
      }
    }
  }
  return 0;
}
