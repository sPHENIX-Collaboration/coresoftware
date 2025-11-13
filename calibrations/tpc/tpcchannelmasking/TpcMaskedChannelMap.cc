
#include "TpcMaskedChannelMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>

#include <ffarawobjects/TpcRawHit.h>
#include <ffarawobjects/TpcRawHitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHPointerListIterator.h>

#include <cdbobjects/CDBTTree.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <trackbase/TpcDefs.h>

#include <TFile.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <regex>


//____________________________________________________________________________..
TpcMaskedChannelMap::TpcMaskedChannelMap(const std::string &name)
  : SubsysReco("TpcMaskedChannelMap")
  , m_run(name)
{
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
int TpcMaskedChannelMap::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator trkr_itr(topNode);
  PHCompositeNode *tpc_node = dynamic_cast<PHCompositeNode *>(trkr_itr.findFirst("PHCompositeNode", "TPC"));
  if (!tpc_node)
  {
    std::cout << __PRETTY_FUNCTION__ << " : ERROR : "
              << "No TPC node found, unregistering module" << std::endl;
    Fun4AllServer::instance()->unregisterSubsystem(this);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHNodeIterator tpc_itr(tpc_node);
  {
    PHPointerListIterator<PHNode> iter(tpc_itr.ls());

    PHNode *thisNode_raw;
    while ((thisNode_raw = iter()))
    {
      if (thisNode_raw->getType() != "PHIODataNode")
      {
        continue;
      }

      PHIODataNode<TpcRawHitContainer> *thisNode = static_cast<PHIODataNode<TpcRawHitContainer> *>(thisNode_raw);
      if (thisNode)
      {
        std::cout << __PRETTY_FUNCTION__ << " : Found TpcRawHitContainer Node "
                  << thisNode->getName() << std::endl;

        TpcRawHitContainer *rawhitcont = (TpcRawHitContainer *) thisNode->getData();
        if (rawhitcont)
        {
          rawhitcont_vec.push_back(rawhitcont);
        }
      }
    }
  }
 
  m_hotFile = "TPCHotMap_" + m_run + ".root";
  m_deadFile = "TPCDeadMap_" + m_run + ".root";

  m_histogramFile = "HIST_" + m_run + ".root";

  h_hits_side0 = new TH2F("h_hits_side0","Hits (Side0);X [cm];Y [cm]",400,-800,800,400,-800,800);
  h_hits_side1 = new TH2F("h_hits_side1","Hits (Side1);X [cm];Y [cm]",400,-800,800,400,-800,800);
  h_masked_side0 = new TH2F("h_masked_side0","Masked Channels (Side0);X [cm];Y [cm]",400,-800,800,400,-800,800);
  h_masked_side1 = new TH2F("h_masked_side1","Masked Channels (Side1);X [cm];Y [cm]",400,-800,800,400,-800,800);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcMaskedChannelMap::process_event(PHCompositeNode * /* unused */)
{
  n_Events += 1.;

  unsigned int raw_hit_num = 0;
  for (TpcRawHitContainer *&rawhitcont : rawhitcont_vec)
  {
    raw_hit_num = rawhitcont->get_nhits();
    for (unsigned int i = 0; i < raw_hit_num; i++)
    {
      auto *hit = rawhitcont->get_hit(i);
      int32_t packet_id = hit->get_packetid();
      int ep = (packet_id - 4000) % 10;
      int sector = (packet_id - 4000 - ep) / 10;
      uint16_t fee = hit->get_fee();
      int channel = hit->get_channel();
      
      if (fee == 16 && channel < 32)
      {
        continue;
      }

      int feeM = FEE_map[fee];
      if (FEE_R[fee] == 2)
      {
        feeM += 6;
      }
      if (FEE_R[fee] == 3)
      {
        feeM += 14;
      }
      double R = M.getR(feeM, channel);
      double phi = 0;
      if (sector < 12)  // NS
      {
        phi = M.getPhi(feeM, channel) + (sector) *M_PI / 6;
      }
      else if (sector >= 12)  // SS
      {
        phi = M.getPhi(feeM, channel) + (sector-12) * M_PI / 6;
      }

      float median = 60;
      for (std::unique_ptr<TpcRawHit::AdcIterator> adc_iterator(hit->CreateAdcIterator());
           !adc_iterator->IsDone();
           adc_iterator->Next())
      {
        const uint16_t sampleN = adc_iterator->CurrentTimeBin();
        const uint16_t adc = adc_iterator->CurrentAdc();
        if (adc - median <= 20)
        {
          continue;
        }

        if (sampleN >= 400 && sampleN <= 430)
        {
          continue;
        }

        nhit_sectors_fees_channels[sector][feeM][channel] += 1;

        if (sector < 12)
        {
          h_hits_side1->Fill(R * cos(phi), R * sin(phi));
        }
        else
        {
          h_hits_side0->Fill(R * cos(phi), R * sin(phi));
        }
      }
    }  
  }
 
  if (raw_hit_num == 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcMaskedChannelMap::End(PHCompositeNode *topNode)
{
  PHG4TpcGeomContainer* geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int nDeadFees = 0;

  int side = 1;
  for (int m_Sector = 0; m_Sector < 24; m_Sector++)
  {
    if (m_Sector > 11)
    {
      side = 0;
    }
    for (int m_FEE = 0; m_FEE < 26; m_FEE++)
    {
      int nDeadInFee = 0;
      for (int m_Channel = 0; m_Channel < 256; m_Channel++)
      {
        //if (nhit_sectors_fees_channels[m_Sector][m_FEE][m_Channel]/n_Events < m_deadChanHitCut)
        if (m_Sector == 1 || m_Sector == 0)
        {
          nDeadInFee++;
          int layer = M.getLayer(m_FEE, m_Channel);

          if (layer < 7)
          {
            continue;
          }

          PHG4TpcGeom* layergeom = geom_container->GetLayerCellGeom(layer);
         
          double phi = 0;
          if (m_Sector < 12)  // NS
          {
            phi = M.getPhi(m_FEE, m_Channel) + (m_Sector) *M_PI / 6;
          }
          else if (m_Sector >= 12)  // SS
          {
            phi = M.getPhi(m_FEE, m_Channel) + (m_Sector-12) * M_PI / 6;
          }
          unsigned int phibin = layergeom->find_phibin(phi,side);
          double R = M.getR(m_FEE, m_Channel);
          m_deadChannelCDB.emplace(layer, mc_sectors[m_Sector % 12], side, phibin,R * cos(phi), R * sin(phi));
          if (side == 0)
          {
            h_masked_side0->Fill(R * cos(phi), R * sin(phi));
          }
          else
          {
            h_masked_side1->Fill(R * cos(phi), R * sin(phi));
          }
        }
        else if (nhit_sectors_fees_channels[m_Sector][m_FEE][m_Channel]/n_Events > m_hotChanHitCut)
        {
          nDeadInFee++;
          int layer = M.getLayer(m_FEE, m_Channel);

          if (layer < 7)
          {
            continue;
          }

          PHG4TpcGeom* layergeom = geom_container->GetLayerCellGeom(layer);
         
          double phi = 0;
          if (m_Sector < 12)  // NS
          {
            phi = M.getPhi(m_FEE, m_Channel) + (m_Sector) *M_PI / 6;
          }
          else if (m_Sector >= 12)  // SS
          {
            phi = M.getPhi(m_FEE, m_Channel) + (m_Sector-12) * M_PI / 6;
          }
          unsigned int phibin = layergeom->find_phibin(phi,side);
          double R = M.getR(m_FEE, m_Channel);
          m_hotChannelCDB.emplace(layer, mc_sectors[m_Sector % 12], side, phibin,R * cos(phi), R * sin(phi));
          if (side == 0)
          {
            h_masked_side0->Fill(R * cos(phi), R * sin(phi));
          }
          else
          {
            h_masked_side1->Fill(R * cos(phi), R * sin(phi));
          }
        }
      }
      if (nDeadInFee > 250)
      {
        nDeadFees++;
      }
    }
  }

  std::cout << "Number of dead FEEs: " << nDeadFees << std::endl;

  CDBTTree cdbttree( m_deadFile );
  cdbttree.SetSingleIntValue( "TotalDeadChannels", m_deadChannelCDB.size() );

  std::cout << "Total Number of Dead Channels: " << m_deadChannelCDB.size() << std::endl;

  int index = 0;
  for( const auto& channel_id:m_deadChannelCDB )
  {
    cdbttree.SetIntValue( index, "layer", channel_id.m_layer );
    cdbttree.SetIntValue( index, "sector", channel_id.m_sector );
    cdbttree.SetIntValue( index, "side", channel_id.m_side );
    cdbttree.SetIntValue( index, "pad", channel_id.m_pad );
    cdbttree.SetFloatValue( index, "x", channel_id.m_x );
    cdbttree.SetFloatValue( index, "y", channel_id.m_y );
    ++index;
  }

  // commit and write
  cdbttree.Commit();
  cdbttree.CommitSingle();
  cdbttree.WriteCDBTTree();  
  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_deadFile << std::endl;

  CDBTTree cdbttree2( m_hotFile );
  cdbttree2.SetSingleIntValue( "TotalHotChannels", m_hotChannelCDB.size() );
  
  std::cout << "Total Number of Hot Channels: " << m_hotChannelCDB.size() << std::endl;

  index = 0;
  for( const auto& channel_id:m_hotChannelCDB )
  {
    cdbttree2.SetIntValue( index, "layer", channel_id.m_layer );
    cdbttree2.SetIntValue( index, "sector", channel_id.m_sector );
    cdbttree2.SetIntValue( index, "side", channel_id.m_side );
    cdbttree2.SetIntValue( index, "pad", channel_id.m_pad );
    cdbttree2.SetFloatValue( index, "x", channel_id.m_x );
    cdbttree2.SetFloatValue( index, "y", channel_id.m_y );
    ++index;
  }

  // commit and write
  cdbttree2.Commit();
  cdbttree2.CommitSingle();
  cdbttree2.WriteCDBTTree();  
  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_hotFile << std::endl;

  TFile *outFile = new TFile(m_histogramFile.c_str(), "RECREATE");
  h_hits_side0->Write();
  h_hits_side1->Write();
  h_masked_side0->Write();
  h_masked_side1->Write();
  outFile->Close();

  delete h_hits_side0; 
  delete h_hits_side1; 
  delete h_masked_side0; 
  delete h_masked_side1;
  delete outFile; 

  return Fun4AllReturnCodes::EVENT_OK;
}
