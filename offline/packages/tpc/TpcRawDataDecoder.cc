
#include "TpcRawDataDecoder.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitSetv1.h>
#include <trackbase/TrkrHitv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TNtuple.h>
#include <limits.h>

#include <string>

//____________________________________________________________________________..
TpcRawDataDecoder::TpcRawDataDecoder(const std::string &name)
  : SubsysReco(name)
 , hm(nullptr)
 , _filename("./outputfile.root")
 {
  std::cout << "TpcRawDataDecoder::TpcRawDataDecoder(const std::string &name)" << std::endl;
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;
  //M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");

  // Open a file, save the ntuple and close the file
  //TFile in_file("/sphenix/user/shulga/Work/Files/pedestal-10616-outfile.root");
  //in_file.GetObject("h_Alive",h_Alive);
  //float chan_id,fee_id,module_id,pedMean,pedStdi, sec_id; 
  //float* row_content;
//
  //  if( Verbosity() )std::cout << "chan_id\t fee_id\t module_id\t pedMean\t pedStdi\t sec_id\n";
  //  for (int irow=0;irow<h_Alive->GetEntries();++irow)
  //  {
  //    h_Alive->GetEntry(irow);
  //    row_content = h_Alive->GetArgs();
  //    chan_id = row_content[0];
  //    fee_id = row_content[1];
  //    module_id = row_content[2];
  //    pedMean = row_content[3];
  //    pedStdi = row_content[4];
  //    sec_id = row_content[5];
  //    if( Verbosity() )
  //    {
  //      std::cout
  //      << chan_id   << "\t"
  //      << fee_id    << "\t" 
  //      << module_id << "\t"
  //      << pedMean   << "\t"
  //      << pedStdi   << "\t"
  //      << sec_id    << "\t"
  //      << std::endl;
  //    }
//
  //    struct ped_tpc_map x
  //    {
  //    };
//
  //    x.CHN_ID = chan_id  ; 
  //    x.FEE_ID = fee_id   ; 
  //    x.MOD_ID = module_id; 
  //    x.PedMean = pedMean  ; 
  //    x.PedStdi = pedStdi  ; 
  //    x.SEC_ID = sec_id   ; 
//
  //    unsigned int key = 256 * (fee_id) + chan_id;
  //    tmap[key] = x;
  //  }
}

//____________________________________________________________________________..
TpcRawDataDecoder::~TpcRawDataDecoder()
{
  delete hm;
  std::cout << "TpcRawDataDecoder::~TpcRawDataDecoder() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int TpcRawDataDecoder::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "TpcRawDataDecoder::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_cdb = CDBInterface::instance();
  std::string calibdir = m_cdb->getUrl("TPC_FEE_CHANNEL_MAP");
  //calibdir = "/sphenix/user/shulga/Work/TpcPadPlane_phi_coresoftware/macros/CDBTest/TPCChannelMap.root";
  if (calibdir[0] == '/')
  {
    // use generic CDBTree to load
    m_cdbttree = new CDBTTree(calibdir.c_str());
    m_cdbttree->LoadCalibrations();
  }
  else
  {
    std::cout << "TpcRawDataDecoder::::InitRun No calibration file found" << std::endl;
    exit(1);
  }
 
  //m_cdbttree = new CDBTTree( "/sphenix/user/shulga/Work/TpcPadPlane_phi_coresoftware/macros/CDBTest/test.root" );
  
  if(m_Debug==1){
    hm = new Fun4AllHistoManager("HITHIST");

    //_h_hit_XYT = new TH3F("_h_hit_XYT" ,"_h_hit_XYT;X, [mm];Y, [mm]; T [BCO]", 400, -800, 800, 400, -800, 800, 500, 128000000000,128050000000);
    _h_hit_XY = new TH2F("_h_hit_XY" ,"_h_hit_XY;X, [mm];Y, [mm]", 400, -800, 800, 400, -800, 800);
    _h_hit_XY_ADCcut = new TH2F("_h_hit_XY_ADCcut" ,"_h_hit_XY_ADCcut;X, [mm];Y, [mm]", 400, -800, 800, 400, -800, 800);
    //_h_hit_PT_ADCcut = new TH3F("_h_hit_PT_ADCcut" ,"_h_hit_PT_ADCcut;Pad number;time [50 ns];BCO;", 400, -0.5, 399.5, 400, -0.5, 399.5, 500, 128000000000,128050000000);

    //hm->registerHisto(_h_hit_XYT );
    //hm->registerHisto(_h_hit_PT_ADCcut);
    hm->registerHisto(_h_hit_XY );
    hm->registerHisto(_h_hit_XY_ADCcut);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawDataDecoder::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TpcRawDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << "TpcRawDataDecoder::InitRun - creating TRKR_HITSET." << std::endl;

    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    hitsetcontainer = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }  
  topNode->print();

  // we reset the BCO for the new run
  starting_BCO = -1;
  rollover_value = 0;
  current_BCOBIN = 0;

  //m_hits = new TrkrHitSetContainerv1();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawDataDecoder::process_event(PHCompositeNode *topNode)
{
  _ievent++;
  //if(_ievent<1270 || _ievent>1270+200) return Fun4AllReturnCodes::DISCARDEVENT;
  // load relevant nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

 PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }


  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  assert( _event );

  if (_event == nullptr)
  {
    std::cout << "TpcRawDataDecoder::Process_Event - Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }
  rollover_value = 0;
  // check all possible TPC packets that we need to analyze
  for(int ep=0;ep<2;ep++){
   for (int sector = 0; sector<24; sector++)
   {
    const int packet = 4000 + sector*10 + ep;

    // Reading packet
    Packet *p = _event->getPacket(packet);

    // Figure out which side
    int side = 0;   
    if(sector>11) side=1;
    if(m_Debug==1){
      char buff[100];
      snprintf(buff, sizeof(buff), "./outputfile_%i.root",sector );
      if(p) _filename = buff;
    }
    if (p)
    {
      std::cout << "My Own TpcRawDataDecoder:: Event getPacket: "<< packet << "| Sector"<< sector << "| EndPoint "<< ep << std::endl;
    }else{
      continue;
    }

    uint64_t triggerBCO = 0;
    int n_tagger = p->lValue(0, "N_TAGGER");
    for (int t = 0; t < n_tagger; t++){
      //  const auto tagger_type = static_cast<uint16_t>(p->lValue(t, "TAGGER_TYPE"));
      // const auto is_endat = static_cast<uint8_t>(p->lValue(t, "IS_ENDAT"));
      const auto is_lvl1 = static_cast<uint8_t>(p->lValue(t, "IS_LEVEL1_TRIGGER"));
      const auto bco = static_cast<uint64_t>(p->lValue(t, "BCO"));
      // const auto lvl1_count = static_cast<uint32_t>(p->lValue(t, "LEVEL1_COUNT"));
      // const auto endat_count = static_cast<uint32_t>(p->lValue(t, "ENDAT_COUNT"));
      // const auto last_bco = static_cast<uint64_t>(p->lValue(t, "LAST_BCO"));
      // const auto modebits = static_cast<uint8_t>(p->lValue(t, "MODEBITS"));
      
      // only printout the is_lvl1 triggers
      // if( Verbosity() )
      if( is_lvl1 )
      {
        triggerBCO = bco;
        std::cout << " is_lvl1: " << (bool) (is_lvl1)
        << " bco: " << bco
        << "_ievent: " << _ievent
        << std::endl;
      }
    }
    //if (triggerBCO != 128330850912) return Fun4AllReturnCodes::DISCARDEVENT;
    int nr_of_waveforms = p->iValue(0, "NR_WF");
    
    //   for (auto &l : m_hitset){
    // l = new TrkrHitSetv1();
    
    int wf;
    //find waveform with earliest BCO in packet
    //take earliest BCO as starting BCO - works only for very long waveform setting
    //current waveforms are ~360 timebins long, eariest BCO == starting BCO should hold for these conditions
    //Check if BCO of previous event is within 360 time bins or 180 BCO. If yes we we are looking at a truncated event and we use the BCO of the previous one as startingBCO

    int earliest_BCO = INT_MAX;

    for (wf = 0; wf < nr_of_waveforms; wf++){      
      int current_BCO = p->iValue(wf, "BCO");
      if(current_BCO < earliest_BCO) earliest_BCO = current_BCO;
      //      std::cout << " earliest BCO:  " << earliest_BCO << " current BCO " << current_BCO << std::endl;
    }
    //  if((earliest_BCO - starting_BCO > 0) && (earliest_BCO - starting_BCO < 180)){
    // starting_BCO = starting_BCO;
    //}else{
    starting_BCO = earliest_BCO;
      //}
    if( Verbosity() )
      std::cout << " _ievent: " << _ievent << " earliest BCO:  " << earliest_BCO << " ep: " << ep << " sector " << sector << " trigBCO: " << triggerBCO << " diff " << triggerBCO - earliest_BCO << std::endl;

    for (wf = 0; wf < nr_of_waveforms; wf++){
      int current_BCO = p->iValue(wf, "BCO") + rollover_value;

      /*
      std::cout << " BCO:  " << p->iValue(wf, "BCO") << std::endl; 
      std::cout << " rollover:  " << rollover_value << std::endl;
      std::cout << " current BCO:  " << current_BCO << std::endl;
      std::cout << " starting BCO:  " << starting_BCO << std::endl;
      std::cout << " earliest BCO:  " << earliest_BCO << std::endl;
      std::cout << " current BCOBIN:  " << current_BCOBIN << std::endl;
      */
      if (starting_BCO < 0)
      {
        starting_BCO = current_BCO;
        }
      
      if (current_BCO < starting_BCO)  // we have a rollover
      {
        rollover_value = 0x100000;
	      current_BCO = p->iValue(wf, "BCO") + rollover_value;
	      starting_BCO = current_BCO;
	      current_BCOBIN++;
      }
      //      std::cout << " nu current BCO:  " << current_BCO << std::endl;
      //std::cout << " nu starting BCO:  " << starting_BCO << std::endl;
      int sampa_nr = p->iValue(wf, "SAMPAADDRESS");
      int channel = p->iValue(wf, "CHANNEL");
      
      int fee = p->iValue(wf, "FEE");
      int samples = p->iValue( wf, "SAMPLES" );
      // clockwise FEE mapping
      //int FEE_map[26]={5, 6, 1, 3, 2, 12, 10, 11, 9, 8, 7, 1, 2, 4, 8, 7, 6, 5, 4, 3, 1, 3, 2, 4, 6, 5};
      int FEE_R[26]={2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
      // conter clockwise FEE mapping (From Takao)
      //int FEE_map[26]={3, 2, 5, 3, 4, 0, 2, 1, 3, 4, 5, 7, 6, 2, 0, 1, 0, 1, 4, 5, 11, 9, 10, 8, 6, 7};
      //int pads_per_sector[3] = {96, 128, 192};
      int FEE_map[26] = { 4,  5,  0,  2,  1,  11,  9,  10,  8,  7,  6,  0,  1,  3,  7,  6,  5,  4,  3,  2,  0,  2,  1,  3,  5,  4};	
      // setting the mapp of the FEE
      int feeM = FEE_map[fee];
      if(FEE_R[fee]==2) feeM += 6;
      if(FEE_R[fee]==3) feeM += 14;
      unsigned int key = 256 * (feeM) + channel;
      //int layer = M.getLayer(feeM, channel);

      std::string varname = "layer";// + std::to_string(key);
      int layer = m_cdbttree->GetIntValue(key,varname);
      // antenna pads will be in 0 layer
      if(layer==0)continue;

      PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
      //varname = "fee" + std::to_string(key);
      //int feeM_cdbt = m_cdbttree->GetIntValue(key,varname);
      
      //varname = "channel" + std::to_string(key);
      //int feeM_chn_cdbt = m_cdbttree->GetIntValue(key,varname);
      
      varname = "R";// + std::to_string(key);
      double R = m_cdbttree->GetDoubleValue(key,varname);
      
      varname = "phi";// + std::to_string(key);
      double phi = pow(-1,side)*m_cdbttree->GetDoubleValue(key,varname) + (sector - side*12 )* M_PI / 6 ;;
      
      //varname = "padN" + std::to_string(key);
      //int padn_cdbt = m_cdbttree->GetIntValue(key,varname);

      //double R = M.getR(feeM, channel);
      //double phi = M.getPhi(feeM, channel) + (sector - side*12 )* M_PI / 6 ;
	    //int pad = M.getPad(feeM, channel);
      //varname = "pad";// + std::to_string(key);
      //int pad = m_cdbttree->GetIntValue(key,varname);

      //std::cout<<"pad "<< padn_cdbt          <<"  "<< pad     << " delta =" << padn_cdbt       - pad     
      //<<"\n  \t  phi "<< phi_cdbt           <<"  "<< phi     << " delta =" << phi_cdbt        - phi     
      //<<"\n  \t  R "<< R_cdbt               <<"  "<< R       << " delta =" << R_cdbt          - R       
      //<<"\n  \t  FEE chn"<< feeM_chn_cdbt   <<"  "<< channel << " delta =" << feeM_chn_cdbt   - channel 
      //<<"\n  \t  FEE "<< feeM_cdbt          <<"  "<< feeM    << " delta =" << feeM_cdbt       - feeM    
      //<<"\n  \t  layer "<< layer_cdbt       <<"  "<< layer   << " delta =" << layer_cdbt      - layer 
      //<< std::endl;        
      int mc_sectors[12] = { 5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
      //int mc_sectors[12] = {5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4};

      float pedestal = 72.4;//round(tmap[key].PedMean);
      TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(layer, (mc_sectors[sector - side*12] ), side);
      TrkrHitSetContainer::Iterator hitsetit = trkrhitsetcontainer->findOrAddHitSet(tpcHitSetKey);
	

      if( Verbosity()>1 ){
	      int sampaAddress = p->iValue(wf, "SAMPAADDRESS");
	      int sampaChannel = p->iValue(wf, "SAMPACHANNEL");
	      int checksum = p->iValue(wf, "CHECKSUM");
	      int checksumError = p->iValue(wf, "CHECKSUMERROR");
	      std::cout << "TpcRawDataDecoder::Process_Event Samples "<< samples 
		        <<" Chn:"<< channel 
		        <<" layer: " << layer 
		        << " sampa: "<< sampa_nr 
		        << " fee: "<< fee 
		        << " Mapped fee: "<< feeM 
		        << " sampaAddress: "<< sampaAddress 
		        << " sampaChannel: "<< sampaChannel 
		        << " checksum: "<< checksum 
		        << " checksumError: "<< checksumError 
		        << " hitsetkey "<< tpcHitSetKey 
		        << " R = " << R
		        << " phi = " << phi
		        << std::endl;

	    }
      pedestal = 0.0;
      for (int s = 0; s < 5; s++)
	    {
	      int adc = p->iValue(wf,s);
	      pedestal += adc;
	    }
      pedestal /=5;
      for (int s = 0; s < samples; s++)
	    {

	      int t = s; // + 2 * (current_BCO - starting_BCO);
	      int adc = p->iValue(wf,s);
	      if(m_Debug==1){
	        _h_hit_XY->Fill(R*cos(phi),R*sin(phi),float(adc)-pedestal);
	      }
	      //	    if(adc-pedestal<4) continue;
	      // generate hit key
	      if(float(adc)-pedestal>2){
          unsigned int phibin = layergeom->find_phibin(phi);
          //if((int)phibin > pads_per_sector[FEE_R[fee]-1]*12) std::cout << "TpcRawDataDecoder:: phibin is out of range > "<< pads_per_sector[FEE_R[fee]-1]*12 << std::endl;
	        TrkrDefs::hitkey hitkey = TpcDefs::genHitKey( phibin, (unsigned int) t);
          //double phi_center = layergeom->get_phicenter(phibin);
          //if(phi_center<0) phi_center += 2*M_PI;
          //int phibin_mc = layergeom->find_phibin(phi);
	        // find existing hit, or create
	        auto hit = hitsetit->second->getHit(hitkey);
  
	        // create hit, assign adc and insert in hitset
	        if (!hit)
		      {
		  
		        // create a new one
		        hit = new TrkrHitv2();
		        hit->setAdc(float(adc)-pedestal);
		        //std::cout<< "ADC = " << adc << " Pedestal = "<< pedestal << "delta = "<< adc-pedestal << std::endl;
		        //if(adc - pedestal > 40) std::cout<< adc - pedestal << "| ";
		        //std::cout<< t << "| ";

		        hitsetit->second->addHitSpecificKey(hitkey, hit);
		      }
	      
	    
	    
	        if(m_Debug==1){
            //if(layer==16 ) _h_hit_PT_ADCcut->Fill(pad, s, triggerBCO,float(adc));
	          if(adc - pedestal > 15){
              //std::cout << "pad = " << pad << "phibin = " << phibin << " phibin_mc = " << phibin_mc << " sector = " << sector << "mc_sectors = " << mc_sectors[sector - side*12] << "phibin/pads_per_sector = " << phibin/pads_per_sector[FEE_R[fee]-1] << " phi_center = " << phi_center << " phi = " << phi << std::endl;
		          _h_hit_XY_ADCcut->Fill(R*cos(phi),R*sin(phi),float(adc)-pedestal);
		          //_h_hit_XYT->Fill(R*cos(phi),R*sin(phi), triggerBCO,float(adc)-pedestal);

	          }
	        }
        }
	      //if (s==samples-1) std::cout << std::endl;
	    }
        
    }
      //     }//auto l:
    
   }// End of run over packets
  }//End of ep loop
  
  // we skip the mapping to real pads at first. We just say
  // that we get 16 rows (segment R2) with 128 pads
  // so each FEE fills 2 rows. Not right, but one step at a time.
  std::cout << "TpcRawDataDecoder:: done" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//int TpcRawDataDecoder::ResetEvent(PHCompositeNode * /*topNode*/)
//{
//  std::cout << "TpcRawDataDecoder::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
//int TpcRawDataDecoder::EndRun(const int runnumber)
//{
//  std::cout << "TpcRawDataDecoder::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
int TpcRawDataDecoder::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "TpcRawDataDecoder::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  if(m_Debug==1)hm->dumpHistos(_filename, "RECREATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
//int TpcRawDataDecoder::Reset(PHCompositeNode * /*topNode*/)
//{
//  std::cout << "TpcRawDataDecoder::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//____________________________________________________________________________..
//void TpcRawDataDecoder::Print(const std::string &what) const
//{
//  std::cout << "TpcRawDataDecoder::Print(const std::string &what) const Printing info for " << what << std::endl;
//}

//____________________________________________________________________________..
//void TpcRawDataDecoder::setHistoFileName(const std::string &what)
//{
//  _filename = what;
//  std::cout << "TpcRawDataDecoder::setHistoFileName(const std::string &what) Histogram File Name is " << what << std::endl;
//}
