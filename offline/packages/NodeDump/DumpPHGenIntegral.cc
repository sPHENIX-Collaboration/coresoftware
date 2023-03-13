#include "DumpPHGenIntegral.h"

#include <phhepmc/PHGenIntegral.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<PHGenIntegral>;

DumpPHGenIntegral::DumpPHGenIntegral(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHGenIntegral::process_Node(PHNode *myNode)
{
  PHGenIntegral *phgenintegral = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phgenintegral = thisNode->getData();
  }
  if (phgenintegral)
  {
    *fout << "get_Integrated_Lumi(): " << phgenintegral->get_Integrated_Lumi() << std::endl;
    *fout << "get_N_Generator_Accepted_Event(): " << phgenintegral->get_N_Generator_Accepted_Event() << std::endl;
    *fout << "get_N_Processed_Event(): " << phgenintegral->get_N_Processed_Event() << std::endl;
    *fout << "get_Sum_Of_Weight(): " << phgenintegral->get_Sum_Of_Weight() << std::endl;
    *fout << "get_CrossSection_Processed_Event(): " << phgenintegral->get_CrossSection_Processed_Event() << std::endl;
    *fout << "get_CrossSection_Generator_Accepted_Event(): " << phgenintegral->get_CrossSection_Generator_Accepted_Event() << std::endl;
    *fout << "get_Description(): " << phgenintegral->get_Description() << std::endl;
  }
  return 0;
}
