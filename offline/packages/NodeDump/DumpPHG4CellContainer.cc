#include "DumpPHG4CellContainer.h"

#include <phool/PHIODataNode.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>

#include <climits>
#include <map>
#include <ostream>
#include <string>
#include <utility>

using MyNode_t = PHIODataNode<PHG4CellContainer>;

DumpPHG4CellContainer::DumpPHG4CellContainer(const std::string &NodeName)
  : DumpObject(NodeName)
{
  return;
}

int DumpPHG4CellContainer::process_Node(PHNode *myNode)
{
  PHG4CellContainer *phg4cellcontainer = nullptr;
  MyNode_t *thisNode = static_cast<MyNode_t *>(myNode);
  if (thisNode)
  {
    phg4cellcontainer = thisNode->getData();
  }
  if (phg4cellcontainer)
  {
    PHG4CellContainer::ConstIterator celler;
    PHG4CellContainer::ConstRange cell_begin_end = phg4cellcontainer->getCells();
    *fout << "size: " << phg4cellcontainer->size() << std::endl;
    for (celler = cell_begin_end.first; celler != cell_begin_end.second; celler++)
    {
      *fout << "id: 0x" << std::hex << celler->second->get_cellid() << std::dec << std::endl;
      if (celler->second->has_binning(PHG4CellDefs::scintillatorslatbinning))
      {
        *fout << "get_column: " << PHG4CellDefs::ScintillatorSlatBinning::get_column(celler->second->get_cellid()) << std::endl;
        *fout << "get_row: " << PHG4CellDefs::ScintillatorSlatBinning::get_row(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::sizebinning))
      {
        *fout << "get_phibin: " << PHG4CellDefs::SizeBinning::get_phibin(celler->second->get_cellid()) << std::endl;
        *fout << "get_zbin: " << PHG4CellDefs::SizeBinning::get_zbin(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::etaphibinning))
      {
        *fout << "get_etabin: " << PHG4CellDefs::EtaPhiBinning::get_etabin(celler->second->get_cellid()) << std::endl;
        *fout << "get_phibin: " << PHG4CellDefs::EtaPhiBinning::get_phibin(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::spacalbinning))
      {
        *fout << "get_etabin: " << PHG4CellDefs::SpacalBinning::get_etabin(celler->second->get_cellid()) << std::endl;
        *fout << "get_phibin: " << PHG4CellDefs::SpacalBinning::get_phibin(celler->second->get_cellid()) << std::endl;
        *fout << "get_fiberid: " << PHG4CellDefs::SpacalBinning::get_fiberid(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::etaxsizebinning))
      {
        *fout << "get_etabin: " << PHG4CellDefs::EtaXsizeBinning::get_etabin(celler->second->get_cellid()) << std::endl;
        *fout << "get_xsizebin: " << PHG4CellDefs::EtaXsizeBinning::get_xsizebin(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::mvtxbinning))
      {
        *fout << "get_index: " << PHG4CellDefs::MVTXBinning::get_index(celler->second->get_cellid()) << std::endl;
      }
      else if (celler->second->has_binning(PHG4CellDefs::tpcbinning))
      {
        *fout << "get_radbin: " << PHG4CellDefs::TPCBinning::get_radbin(celler->second->get_cellid()) << std::endl;
        *fout << "get_phibin: " << PHG4CellDefs::TPCBinning::get_phibin(celler->second->get_cellid()) << std::endl;
      }
      else
      {
        *fout << "binning "
              << PHG4CellDefs::get_binning(celler->second->get_cellid())
              << " for detid: "
              << PHG4CellDefs::get_detid(celler->second->get_cellid())
              << " not implemented in DumpPHG4CellContainer" << std::endl;
      }

      for (auto ic = 0; ic < UCHAR_MAX; ic++)
      {
        PHG4Cell::PROPERTY prop_id = static_cast<PHG4Cell::PROPERTY>(ic);
        if (celler->second->has_property(prop_id))
        {
          *fout << "prop id: " << static_cast<unsigned int>(ic);
          std::pair<const std::string, PHG4Cell::PROPERTY_TYPE> property_info = PHG4Cell::get_property_info(prop_id);
          *fout << ", name " << property_info.first << " value ";
          switch (property_info.second)
          {
          case PHG4Cell::type_int:
            *fout << celler->second->get_property_int(prop_id);
            break;
          case PHG4Cell::type_uint:
            *fout << celler->second->get_property_uint(prop_id);
            break;
          case PHG4Cell::type_float:
            *fout << celler->second->get_property_float(prop_id);
            break;
          default:
            *fout << " unknown type ";
          }
          *fout << std::endl;
        }
      }
      PHG4Cell::EdepConstRange hitedep_begin_end = celler->second->get_g4hits();
      for (PHG4Cell::EdepConstIterator iter = hitedep_begin_end.first; iter != hitedep_begin_end.second; ++iter)
      {
        *fout << "hit 0x" << std::hex << iter->first << std::dec << " edep: " << iter->second << std::endl;
      }
      PHG4Cell::ShowerEdepConstRange shower_begin_end = celler->second->get_g4showers();
      for (PHG4Cell::ShowerEdepConstIterator iter = shower_begin_end.first; iter != shower_begin_end.second; ++iter)
      {
        *fout << "shower 0x" << std::hex << iter->first << std::dec << " edep: " << iter->second << std::endl;
      }
    }
  }
  return 0;
}
