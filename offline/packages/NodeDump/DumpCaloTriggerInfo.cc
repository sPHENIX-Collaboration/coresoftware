#include "DumpCaloTriggerInfo.h"

#include <calotrigger/CaloTriggerInfo.h>

#include <phool/PHIODataNode.h>

#include <ostream>
#include <string>

using MyNode_t = PHIODataNode<CaloTriggerInfo>;

DumpCaloTriggerInfo::DumpCaloTriggerInfo(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpCaloTriggerInfo::process_Node(PHNode *myNode)
{
  CaloTriggerInfo *calotriggerinfo = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    calotriggerinfo = thisNode->getData();
  }
  if (calotriggerinfo)
  {
    *fout << " get_best_EMCal_2x2_E: " << calotriggerinfo->get_best_EMCal_2x2_E() << std::endl;
    *fout << " get_best_EMCal_2x2_eta: " << calotriggerinfo->get_best_EMCal_2x2_eta() << std::endl;
    *fout << " get_best_EMCal_2x2_phi: " << calotriggerinfo->get_best_EMCal_2x2_phi() << std::endl;
    *fout << " get_best_EMCal_4x4_E: " << calotriggerinfo->get_best_EMCal_4x4_E() << std::endl;
    *fout << " get_best_EMCal_4x4_eta: " << calotriggerinfo->get_best_EMCal_4x4_eta() << std::endl;
    *fout << " get_best_EMCal_4x4_phi: " << calotriggerinfo->get_best_EMCal_4x4_phi() << std::endl;
    *fout << " get_best2_EMCal_4x4_E: " << calotriggerinfo->get_best2_EMCal_4x4_E() << std::endl;
    *fout << " get_best2_EMCal_4x4_eta: " << calotriggerinfo->get_best2_EMCal_4x4_eta() << std::endl;
    *fout << " get_best2_EMCal_4x4_phi: " << calotriggerinfo->get_best2_EMCal_4x4_phi() << std::endl;
    *fout << " get_best_FullCalo_0p2x0p2_E: " << calotriggerinfo->get_best_FullCalo_0p2x0p2_E() << std::endl;
    *fout << " get_best_FullCalo_0p2x0p2_eta: " << calotriggerinfo->get_best_FullCalo_0p2x0p2_eta() << std::endl;
    *fout << " get_best_FullCalo_0p2x0p2_phi: " << calotriggerinfo->get_best_FullCalo_0p2x0p2_phi() << std::endl;
    *fout << " get_best_FullCalo_0p4x0p4_E: " << calotriggerinfo->get_best_FullCalo_0p4x0p4_E() << std::endl;
    *fout << " get_best_FullCalo_0p4x0p4_eta: " << calotriggerinfo->get_best_FullCalo_0p4x0p4_eta() << std::endl;
    *fout << " get_best_FullCalo_0p4x0p4_phi: " << calotriggerinfo->get_best_FullCalo_0p4x0p4_phi() << std::endl;
    *fout << " get_best_FullCalo_0p6x0p6_E: " << calotriggerinfo->get_best_FullCalo_0p6x0p6_E() << std::endl;
    *fout << " get_best_FullCalo_0p6x0p6_eta: " << calotriggerinfo->get_best_FullCalo_0p6x0p6_eta() << std::endl;
    *fout << " get_best_FullCalo_0p6x0p6_phi: " << calotriggerinfo->get_best_FullCalo_0p6x0p6_phi() << std::endl;
    *fout << " get_best_FullCalo_0p8x0p8_E: " << calotriggerinfo->get_best_FullCalo_0p8x0p8_E() << std::endl;
    *fout << " get_best_FullCalo_0p8x0p8_eta: " << calotriggerinfo->get_best_FullCalo_0p8x0p8_eta() << std::endl;
    *fout << " get_best_FullCalo_0p8x0p8_phi: " << calotriggerinfo->get_best_FullCalo_0p8x0p8_phi() << std::endl;
    *fout << " get_best_FullCalo_1p0x1p0_E: " << calotriggerinfo->get_best_FullCalo_1p0x1p0_E() << std::endl;
    *fout << " get_best_FullCalo_1p0x1p0_eta: " << calotriggerinfo->get_best_FullCalo_1p0x1p0_eta() << std::endl;
    *fout << " get_best_FullCalo_1p0x1p0_phi: " << calotriggerinfo->get_best_FullCalo_1p0x1p0_phi() << std::endl;
  }
  return 0;
}
