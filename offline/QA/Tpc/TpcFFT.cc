// Includes
#include "TpcFFT.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <cassert>
#include <cstddef>
#include <memory>

#include <qautils/QAHistManagerDef.h>
#include <boost/format.hpp>

#include <iostream>
#include <string>

#include </sphenix/user/llegnosky/TPCAnalysis/TpcModules/Biquad.h>
//

//____________________________________________________________________________..
TpcFFT::TpcFFT(const std::string &name)
  : SubsysReco(name)
{
// reserves memory for max ADC samples
m_adcSamples.resize(1024, 0);
}
  
//____________________________________________________________________________..
int TpcFFT::InitRun(PHCompositeNode * /*unused*/)
{

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcFFT::process_event(PHCompositeNode *topNode)
{
  char name[100];
  
  // Defines object from class Event which calls getClass function from
  // findNode class
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");

  // Checks if event exists and returns error if not
  if (_event == nullptr)
    {
      std::cout << "TPCRawDataTree::Process_Event - Event not found" << std::endl;
      return -1;
    }

  // Checks if event is "special" and discards it if so
  if (_event->getEvtType() >= 8)  /// special events
    {
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
//

  // Initialize HistoManager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // Reference histograms initialized in header file to histos in HistoManager
  h_WF = dynamic_cast<TH2F*>(hm->getHisto(boost::str(boost::format("h_WF_Dist")).c_str()));
  h_FFT = dynamic_cast<TH2F*>(hm->getHisto(boost::str(boost::format("h_FFT_Dist")).c_str()));

  Biquad *lpFilter = new Biquad();

  lpFilter->setBiquad(bq_type_lowpass, 9.0/20.0, 0.707, 0);

  std::vector<Packet *> pktvec = _event->getPacketVector();

// Loop over packets in event
for (auto packet : pktvec)
  {

    int32_t packet_id = packet->getIdentifier();

    if (Verbosity())
      {
	std::cout << __PRETTY_FUNCTION__ << " : decoding packet " << packet_id << std::endl;
      }

    if (!packet)
    {
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " : missing packet " << packet_id << std::endl;
      }
      continue;
    }

    int ep = (packet_id-4000) % 10;
    sector = (packet_id - 4000 - ep)/10;

    // pull number of waveforms
    m_nWaveformInFrame = packet->iValue(0, "NR_WF");

    for (int wf = 0; wf < m_nWaveformInFrame; wf++)
      {
	
	// // Create and register histos in HistoManager
	// {
	//   auto h = new TH1F(boost::str(boost::format("h_WF_%s_%d") % sec[sector].c_str() % ent_num).c_str(),";t bins [50 ns];ADC",360,0,360);
	//   hm->registerHisto(h);
	// }

	// {
	//   auto h = new TH1F(boost::str(boost::format("h_FFT_%s_%d") % sec[sector].c_str() % ent_num).c_str(),";Frequency [20*bin/360];Counts",360,0,20);
	//   hm->registerHisto(h);
	// }

	// // Reference histograms initialized in header file to histos in HistoManager
	// h_WF = dynamic_cast<TH1F*>(hm->getHisto(boost::str(boost::format("h_WF_%s_%d") % sec[sector].c_str() % ent_num).c_str()));
	// h_FFT = dynamic_cast<TH1F*>(hm->getHisto(boost::str(boost::format("h_FFT_%s_%d") % sec[sector].c_str() % ent_num).c_str()));

	// Checks if sample number and number of ADC values agrees
	//assert(m_nSamples < (int) m_adcSamples.size());
	if(m_nSamples > (int) m_adcSamples.size()){
	  //ent_num++;
	  continue;
	}
	

	dead = false;
	// Loop over samples in waveform
	for (int s = 0; s < m_nSamples; s++)
	  {
	    // Assign ADC value of sample s in waveform wf to adcSamples[s]
	    m_adcSamples[s] = packet->iValue(wf, s);
	
	    if(m_adcSamples[s] == 0 || TMath::IsNaN(float(m_adcSamples[s]))) {
	      dead=true;
	      break;
	    }
	    if(s<360) Samples[s]=m_adcSamples[s];
	  }

	if(dead){
	  //ent_num++;
	  continue;
	}

	sprintf(name,"h_WF_%s_%d",sec[sector].c_str(),ent_num);
	TH1F *h_WF_temp = new TH1F(name,name,360,0,360);
	
	sprintf(name,"h_FFT_%s_%d",sec[sector].c_str(),ent_num);
	TH1 *h_FFT_temp = new TH1F(name,name,360,0,20);
	
	m_FEE = packet->iValue(wf, "FEE");
	m_Channel = packet->iValue(wf, "CHANNEL");
	m_nSamples = packet->iValue(wf, "SAMPLES");

	pedestal.push_back(TMath::Median(360,Samples));

	for(int adc_sam_no=0;adc_sam_no<m_nSamples;adc_sam_no++){
	  if(m_adcSamples[adc_sam_no]<1024){
	    h_WF_temp->Fill(adc_sam_no,m_adcSamples[adc_sam_no]);
	    if(m_adcSamples[adc_sam_no]>pedestal[wf]+50.0 || m_adcSamples[adc_sam_no]<pedestal[wf]-50.0) continue;
	    pedestal_sum += m_adcSamples[adc_sam_no];
	    pedestal_sigma_sum += pow(m_adcSamples[adc_sam_no],2);
	    pedestal_samples += 1.0;
	  }
	}

	pedestal_sigma.push_back(sqrt(pedestal_sigma_sum/pedestal_samples - pow(pedestal_sum/pedestal_samples,2)));

	h_FFT_temp = h_WF_temp->FFT(h_FFT_temp,"MAG");
	Int_t FFT_entries = h_FFT_temp->GetEntries();
	h_FFT_temp->GetXaxis()->SetRangeUser(20.0/360.0,10.0);
	//h_FFT_temp->Scale(1.0/sqrt(FFT_entries));

	Int_t max_f_bin = h_FFT_temp->GetMaximumBin();
	Float_t max_f = (max_f_bin-1)*20.0/360.0;

	// if(max_f>8.0){
	//   for(int ent=0;ent<FFT_entries;ent++){
	//     Float_t old_bin_value = h_FFT->GetBinContent(ent);
	//     Float_t new_bin_value = lpFilter->process(old_bin_value);
	//     h_FFT->SetBinContent(ent,new_bin_value);
	//   }
	// }

	if(max_f>8.0){
	  for(int ent=0;ent<FFT_entries;ent++){
	    Float_t old_bin_value = h_FFT_temp->GetBinContent(ent);
	    Float_t new_bin_value = lpFilter->process(old_bin_value);
	    h_FFT_temp->SetBinContent(ent,new_bin_value);
	  }
	  for(int sam=0;sam<360;sam++){
	    h_WF->Fill(sam,h_WF_temp->GetBinContent(sam+1));
	    h_FFT->Fill(sam*20.0/360.0,h_FFT_temp->GetBinContent(sam+1)*1.0/sqrt(FFT_entries));
	  }
	}

	// if(max_f>8.0){
	//   TH1 *WF_temp;
	//   WF_temp=(TH1*)h_WF->Clone();
	//   WF_temp->SetDirectory(0);
	//   WF_clone.push_back(WF_temp);
	//   TH1 *FFT_temp;
	//   FFT_temp=(TH1*)h_FFT->Clone();
	//   FFT_temp->SetDirectory(0);
	//   FFT_clone.push_back(FFT_temp);
	//   evt_num.push_back(ent_num);
	// }
	
	ent_num++;
      }
  }

 // sprintf(name, "/sphenix/user/llegnosky/TPCAnalysis/RootFiles/TpcFFT_QA_29746.root");
 // std::string outfile = name;
 // TpcFFTfile = new TFile(outfile.c_str(),"RECREATE");

 // int length = FFT_clone.size();
 
 // sprintf(name,"Event_Numbers_sec%s",sec[sector].c_str());
 // TpcFFTfile->WriteObject(&evt_num,name);
 // for(int evt=0;evt<length;evt++){
 //   sprintf(name,"WF_sec%s_evt%d",sec[sector].c_str(),evt_num[evt]);
 //   WF_clone[evt]->Write(name);
 //   sprintf(name,"FFT_sec%s_evt%d",sec[sector].c_str(),evt_num[evt]);
 //   FFT_clone[evt]->Write(name);
 //   sprintf(name,"Ped_sec%s_evt%d",sec[sector].c_str(),evt_num[evt]);
 //   TpcFFTfile->WriteObject(&pedestal,name);
 //   sprintf(name,"Ped_RMS_sec%s_evt%d",sec[sector].c_str(),evt_num[evt]);
 //    TpcFFTfile->WriteObject(&pedestal_sigma,name);
 // }
 
 // TpcFFTfile->Close();
 
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcFFT::End(PHCompositeNode * /*unused*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
void TpcFFT::createHistos()
{
  // Initialize HistoManager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // Create and register histos in HistoManager
  {
    auto h = new TH2F(boost::str(boost::format("h_WF_Dist")).c_str(),";t bins [50 ns];ADC",360,0,360,1024,0,1024);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(boost::str(boost::format("h_FFT_Dist")).c_str(),";Frequency [20*bin/360];Counts",360,0,20,500,0,1000);
    hm->registerHisto(h);
  }
}
