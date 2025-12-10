#include "DiodeReco.h"

#include <ffamodules/CDBInterface.h>

#include <ffarawobjects/TpcDiodeContainerv1.h>
#include <ffarawobjects/TpcDiodev1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/caen_correction.h>
#include <Event/packet.h>


#include <TSystem.h>

#include <algorithm>
#include <cstdint>                             // for uint16_t
#include <cstdlib>                             // for exit, size_t
#include <iostream>                             // for basic_ostream, operat...

DiodeReco::DiodeReco(const std::string &name)
  : SubsysReco(name)
  , m_DiodeContainerName("TPCDIODES")
{
  for (int c = 0; c < 32; c++)
  {
    adc.clear();
  }
}

int DiodeReco::InitRun(PHCompositeNode *topNode)
{
   std::cout << "DiodeReco::InitRun(PHCompositeNode *topNode) Initializing" << std::endl;

  m_cdb = CDBInterface::instance();
  
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "TPC"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("TPC");
    dstNode->addNode(detNode);
  }

  TpcDiodeContainer *tpcdiodecont = findNode::getClass<TpcDiodeContainer>(detNode, m_DiodeContainerName);
  if (!tpcdiodecont)
  {
    tpcdiodecont = new TpcDiodeContainerv1();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(tpcdiodecont, m_DiodeContainerName, "PHObject");
    detNode->addNode(newNode);
  }

  // char name[100];
  // char title[100];

  // c_persistency_N = new TCanvas("c_persistency_N","c_persistency_N",800,600);
  // c_persistency_N->Divide(4,4);
  // c_persistency_S = new TCanvas("c_persistency_S","c_persistency_S",800,600);
  // c_persistency_S->Divide(4,4);
  // c_waveforms = new TCanvas("c_waveforms","c_waveforms",800,600);

  // for(int c=0;c<32;c++){
  //   sprintf(name,"CAEN Persistency Waveform (Channel %d)",c);
  //   sprintf(title,"CAEN Persistency Waveform (Channel %d);time bins [ns];ADC",c);
  //   persistency[c] = new TH2F(name,title,1024,-0.5,1023.5,256,-250,4250);
  // }
  // waveforms = new TH2F("CAEN Waveforms","CAEN Waveforms;time bins [ns];channels",1024,0,1024,32,0,32);
  return Fun4AllReturnCodes::EVENT_OK;
}

int DiodeReco::process_event(PHCompositeNode *topNode)
{
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (!_event)
  {
    std::cout << "RawdataUnpacker::Process_Event - Event not found" << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  diodes = findNode::getClass<TpcDiodeContainer>(topNode, "TPCDIODES");

  //_event->identify();
  caen_correction *cc{nullptr};
  Packet *p = _event->getPacket(2000);
  if (p)
  {
    if (!cc)
    {
      switch (p->iValue(0, "FREQUENCY"))
      {
      case 0:  // 5GS
	calibdir = m_cdb->getUrl("TPC_CAEN_CORRECTION_24974_5G");
        cc = new caen_correction(calibdir.c_str());
        break;
      case 1:  // 2.5GS
        calibdir = m_cdb->getUrl("TPC_CAEN_CORRECTION_24974_2G");
        cc = new caen_correction(calibdir.c_str());
        break;
      case 2:  // 1GS
        calibdir = m_cdb->getUrl("TPC_CAEN_CORRECTION_24974_1G");
        cc = new caen_correction(calibdir.c_str());
        break;
      default:
        std::cout << "Bad selection " << p->iValue(0, "FREQUENCY") << std::endl;
        gSystem->Exit(1);
        exit(1);
        break;
      }
    }

    if (cc)
    {
      cc->init(p);
    }

    TpcDiodev1 *newdiode = new TpcDiodev1();

    for (int c = 0; c < 32; c++)
    {
      for (int t = 0; t < 1024; t++)
      {
        adc.push_back(cc->caen_corrected(t, c));
      }

      PedestalCorrected(0, 200);
      uint16_t maxadc = MaxAdc(2, 200, 300);
      uint16_t maxbin = MaxBin(2);
      double integral = Integral(0, 1024);
      uint16_t nabovethreshold = NAboveThreshold(200, 100);
      double pulsewidth = PulseWidth(200, 100);
      const uint16_t samples = adc.size();

      uint16_t packetid = 2000;
      uint16_t channel = c;

      newdiode->set_packetid(packetid);
      newdiode->set_channel(channel);
      newdiode->set_maxadc(maxadc);
      newdiode->set_maxbin(maxbin);
      newdiode->set_integral(integral);
      newdiode->set_nabovethreshold(nabovethreshold);
      newdiode->set_pulsewidth(pulsewidth);
      newdiode->set_samples(samples);

      for (uint16_t s = 0; s < samples; s++)
      {
        uint16_t adcval = adc[s];
        newdiode->set_adc(s, adcval);
      }

      diodes->AddDiode(newdiode);
      adc.clear();
    }

    // if(event==2){
    // 	  for(int c=0;c<32;c++)
    // 	    {
    // 	      for(int t=0;t<1024;t++){
    // 		adc.push_back(cc->caen_corrected(t,c));
    // 	      }

    // 	      PedestalCorrected(0,200);

    // 	      if(nlaser>1)
    // 		{
    // 		  std::cout << "More than one laser fired! Stopping code!" << std::endl;
    // 		  return -1;
    // 		}

    // 	      if(laser<0)
    // 		{
    // 		  std::cout << "No laser fired in this event!" << std::endl;
    // 		}
    // 	      for(int s=0;s<static_cast<int>(adc.size());s++){
    // 		persistency[c]->Fill(s,adc[s]);
    // 		waveforms->Fill(s,c,adc[s]);
    // 	      }
    // 	      if(c<16){
    // 		c_persistency_N->cd(c+1);
    // 		persistency[c]->Draw("colz");
    // 	      }
    // 	      else
    // 		{
    // 		  c_persistency_S->cd(c-15);
    // 		  persistency[c]->Draw("colz");
    // 		}
    // 	      adc.clear();
    // 	    }
    // }
    // p->dump();
    std::cout << diodes->get_Laser() << std::endl;
    // for(int c=0;c<32;c++){
    // 	std::cout << diodes->get_diode(c)->get_maxadc() << std::endl;
    // }
    delete p;
  }

  // c_waveforms->cd();
  // waveforms->Draw("LEGO");

  return Fun4AllReturnCodes::EVENT_OK;
}

double DiodeReco::MaxAdc(int n, int low_bin, int high_bin) const
{
  double MaxSum = -99999;  // Maximum sum over n bins within the bin range
  int MaxMid = -1;         // Bin number of the middle bin used to calculate MaxSum

  // Bracket the limits safely...
  int start = std::max(1 + n, low_bin + n);
  int end = std::min(int(adc.size()) - n, high_bin - n);
  for (int mid = start; mid < end; mid++)
  {
    double Sum = 0;
    for (int bin = mid - n; bin <= mid + n; bin++)
    {
      Sum += adc[bin];
    }
    if (Sum > MaxSum)
    {
      MaxSum = Sum;
      MaxMid = mid;
    }
  }

  if (MaxMid < 0)
  {
    std::cout << "Error: Maximum ADC could not be found!" << std::endl;
  }
  return MaxSum / (2.0 * n + 1.0);
}

int DiodeReco::MaxBin(int n) const
{
  double MaxSum = -99999;
  int MaxMid = -1;
  for (int mid = 1 + n; mid < static_cast<int>(adc.size()) - n; mid++)
  {
    double Sum = 0;
    for (int bin = mid - n; bin <= mid + n; bin++)
    {
      Sum += adc[bin];
    }
    if (Sum > MaxSum)
    {
      MaxSum = Sum;
      MaxMid = mid;
    }
  }
  if (MaxMid < 0)
  {
    std::cout << "Error: Maximum ADC could not be found!" << std::endl;
  }
  return MaxMid;
}

double DiodeReco::Integral(int low_bin, int high_bin) const
{
  low_bin = std::max(low_bin, 0);
  high_bin = std::min(high_bin, static_cast<int>(adc.size()));

  double SUM = 0;
  for (int i = low_bin; i < high_bin; i++)
  {
    SUM += adc[i];
  }
  return SUM;
}

int DiodeReco::NAboveThreshold(double upper_thr, double lower_thr) const
{
  int nAbove = 0;

  bool belowThreshold = true;

  for (double adc_val : adc)
  {
    if (belowThreshold && adc_val >= upper_thr)
    {
      nAbove++;
      belowThreshold = false;
    }

    else if (!belowThreshold && adc_val < lower_thr)
    {
      belowThreshold = true;
    }
  }

  return nAbove;
}

double DiodeReco::PulseWidth(double upper_thr, double lower_thr) const
{
  //  The results of this routine are ONLY valid
  //  if NAbove is one.

  bool belowThreshold = true;

  int left = 0;
  int right = 0;

  for (int i = 0; i < static_cast<int>(adc.size()); i++)
  {
    if (belowThreshold && adc[i] >= upper_thr)
    {
      left = i;
      belowThreshold = false;
    }

    else if (!belowThreshold && adc[i] < lower_thr)
    {
      right = i;
      belowThreshold = true;
    }
  }

  return right - left;
}

void DiodeReco::PedestalCorrected(int low_bin = -1, int high_bin = 9999)
{
  low_bin = std::max(low_bin, 0);
  if (high_bin > static_cast<int>(adc.size()))
  {
    high_bin = adc.size();
  }

  // Copy all voltages in the pedestal region & sort
  std::vector<double> sam;
  for (int i = low_bin; i < high_bin; i++)
  {
    sam.push_back(adc[i]);
  }
  sort(sam.begin(), sam.end());

  // Assign the pedestal as the median of this distribution
  int n = sam.size();
  double PEDESTAL;
  if (n % 2 != 0)
  {
    PEDESTAL = sam[n / 2];
  }
  else
  {
    PEDESTAL = (sam[(n - 1) / 2] + sam[n / 2]) / 2.0;
  }

  for (double & adc_val : adc)
  {
    adc_val = adc_val - PEDESTAL;
  }
}
