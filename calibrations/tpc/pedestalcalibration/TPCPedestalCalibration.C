#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

/*************************************************************/
/*                TPC Pedestal Calibration                   */
/*               Thomas Marshall,Aditya Dash                 */
/*        rosstom@ucla.edu,aditya55@physics.ucla.edu         */
/*************************************************************/

void TPCPedestalCalibration(vector<string> inputFiles = {"/sphenix/user/rosstom/test/testFiles/TPC_ebdc00_pedestal-00010131-0000.prdf_TPCRawDataTree.root"}){
  for (int fileNum = 0; fileNum < inputFiles.size(); fileNum++){
    string sectorNum = inputFiles[fileNum];
    size_t pos = sectorNum.find("ebdc");
    sectorNum.erase(pos+6);
    sectorNum.erase(sectorNum.begin(),sectorNum.begin()+pos+4); 
    Int_t sector;
    if (sectorNum.substr(0,1) == "0")
    {
      sectorNum.erase(sectorNum.begin(),sectorNum.end()-1);
      sector = std::stoi(sectorNum);
    }
    else sector = std::stoi(sectorNum); 
 
    Int_t mod_arr[26]={2,2,1,1,1,3,3,3,3,3,3,2,2,1,2,2,1,1,2,2,3,3,3,3,3,3};
    Int_t slot_arr[26] = {5,6,1,3,2,12,10,11,9,8,7,1,2,4,8,7,6,5,4,3,1,3,2,4,6,5}; 

    TFile *pedestalInput = TFile::Open(inputFiles[fileNum].c_str());
    TTreeReader myReader("SampleTree", pedestalInput);

    TTreeReaderValue<Int_t> nSamples_in(myReader, "nSamples");
    TTreeReaderValue<Int_t> fee_in(myReader, "fee");
    TTreeReaderArray<UShort_t> adcSamples_in(myReader, "adcSamples");
    TTreeReaderValue<Int_t> Channel_in(myReader, "Channel");

//Adding the output tree
    TString* outputfilename=new TString("pedestalCalibration_TPC_ebdc"+sectorNum+".root");
    TFile* outputfile=new TFile(outputfilename->Data(),"recreate");
    TTree* outputTree=new TTree("outputTree","outputTree");
    Int_t isAlive=1; // 1 if the channel is working properly, 0 if no signal(number of samples is 0, adc value is 0 or nan), pedestal above 200 or below 10
    Float_t pedMean; // average pedestal value over all samples
    Float_t pedStd; // pedestal standard deviation over all samples
    Int_t chan; // channel number
    Int_t fee; // fee number
    Int_t module; // module number (e.g. R1, R2, R3)
    Int_t slot; // slot number
    outputTree->Branch("isAlive",&isAlive,"isAlive/I"); 
    outputTree->Branch("pedMean",&pedMean,"pedMean/F");
    outputTree->Branch("pedStd",&pedStd,"pedStd/F");
    outputTree->Branch("sector",&sector,"sector/I");
    outputTree->Branch("fee",&fee,"fee/I");
    outputTree->Branch("chan",&chan,"chan/I");
    outputTree->Branch("module",&module,"module/I");
    outputTree->Branch("slot",&slot,"slot/I");

    Float_t ave_adc_fee_channel[26][256];
    Float_t std_adc_fee_channel[26][256];
    Float_t counts_adc_fee_channel[26][256];
    Int_t alive_array_fee_channel[26][256];
    for(int fee_no=0;fee_no<26;fee_no++){
        for(int channel_no=0;channel_no<256;channel_no++)
	{
        	ave_adc_fee_channel[fee_no][channel_no]=0.0;
        	std_adc_fee_channel[fee_no][channel_no]=0.0;
        	counts_adc_fee_channel[fee_no][channel_no]=0.0;
		alive_array_fee_channel[fee_no][channel_no]=1;
        }
    }

    while(myReader.Next())
        {
            if(*nSamples_in==0) 
	    {
		alive_array_fee_channel[*fee_in][*Channel_in]=0; 
		continue;
	    }
            
            bool dead = false;
            for(int adc_sam_no=0;adc_sam_no<*nSamples_in;adc_sam_no++)
            {
                if (adcSamples_in[adc_sam_no] == 0 || TMath::IsNaN(float(adcSamples_in[adc_sam_no]))) 
                {
			alive_array_fee_channel[*fee_in][*Channel_in]=0;
		}
            }
            if (dead) {continue;}

            for(int adc_sam_no=0;adc_sam_no<*nSamples_in;adc_sam_no++){
                ave_adc_fee_channel[*fee_in][*Channel_in]+=adcSamples_in[adc_sam_no];
                std_adc_fee_channel[*fee_in][*Channel_in]+=pow(adcSamples_in[adc_sam_no],2);
                counts_adc_fee_channel[*fee_in][*Channel_in]+=1.0;
            }

        }

    for(int fee_no=0;fee_no<26;fee_no++)
    {
          for(int channel_no=0;channel_no<256;channel_no++)
          {
            if(counts_adc_fee_channel[fee_no][channel_no] != 0.0) 
            {
                float temp1 = ave_adc_fee_channel[fee_no][channel_no]/counts_adc_fee_channel[fee_no][channel_no];
                float temp2 = std_adc_fee_channel[fee_no][channel_no]/counts_adc_fee_channel[fee_no][channel_no];
                ave_adc_fee_channel[fee_no][channel_no] = temp1;
                std_adc_fee_channel[fee_no][channel_no] = temp2;
            }
            else
	    {
                ave_adc_fee_channel[fee_no][channel_no] = 0.0;
                std_adc_fee_channel[fee_no][channel_no] = 0.0;
		alive_array_fee_channel[fee_no][channel_no]=0;
            }
        
            if(ave_adc_fee_channel[fee_no][channel_no] > 200 || ave_adc_fee_channel[fee_no][channel_no] < 10)
            {
		alive_array_fee_channel[fee_no][channel_no]=0;
            }

	    pedMean=ave_adc_fee_channel[fee_no][channel_no];
	    pedStd=sqrt(std_adc_fee_channel[fee_no][channel_no] - pow(ave_adc_fee_channel[fee_no][channel_no],2));
	    isAlive=alive_array_fee_channel[fee_no][channel_no];
	    chan=channel_no;
            fee=fee_no;
	    module=mod_arr[fee_no];
	    slot=slot_arr[fee_no];
	    outputTree->Fill(); 
          }
      }
      outputfile->cd();
      outputTree->Write();
      outputfile->Close();
      pedestalInput->Close(); 
  }
}
