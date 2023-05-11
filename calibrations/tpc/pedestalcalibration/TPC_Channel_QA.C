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


void TPC_Channel_QA(){
  gROOT->SetBatch(kTRUE);
  std::ofstream outdata;
  outdata.open("noisyChannels.txt");
  if( !outdata ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
  }
  for (int q = 0; q < 24; q++){
    string sectorNumber;
    if (q < 10) sectorNumber = "0"+std::to_string(q);
    else sectorNumber = std::to_string(q);   

    string runNumber = "pedestal-00010131";
    
    string fileName = "/sphenix/user/rosstom/test/testFiles/TPC_ebdc"+sectorNumber+"_"+runNumber+"-0000.prdf_TPCRawDataTree.root";
 
//    string mod[26] = {"R2","R2","R1","R1","R1","R3","R3","R3","R3","R3","R3","R2","R2","R1","R2","R2","R1","R1","R2","R2","R3","R3","R3","R3","R3","R3"}; 
//    string slot[26] = {"5","6","1","3","2","12","10","11","9","8","7","1","2","4","8","7","6","5","4","3","1","3","2","4","6","5"}; 
    int mod_arr[26]={2,2,1,1,1,3,3,3,3,3,3,2,2,1,2,2,1,1,2,2,3,3,3,3,3,3};
    int slot_arr[26] = {5,6,1,3,2,12,10,11,9,8,7,1,2,4,8,7,6,5,4,3,1,3,2,4,6,5}; 

    TFile *truthInput = TFile::Open(fileName.c_str());
    TTreeReader myReader("SampleTree", truthInput);

    TTreeReaderValue<Int_t> nSamples_in(myReader, "nSamples");
    TTreeReaderValue<Int_t> fee_in(myReader, "fee");
    TTreeReaderArray<UShort_t> adcSamples_in(myReader, "adcSamples");
    TTreeReaderValue<Int_t> Channel_in(myReader, "Channel");

//Adding the output tree
        TString* outputfilename=new TString("outputfile_TPC_ebdc"+sectorNumber+"_"+runNumber+".root");
 	TFile* outputfile=new TFile(outputfilename->Data(),"recreate");
 	TTree* outputTree=new TTree("outputTree","outputTree");
	Int_t isAlive=1.0; // 1 if the channel is working properly, 0 if no signal(number of samples is 0, adc value is 0 or nan), pedestal above 200 or below 10
	Float_t pedMean; // average pedestal value over all samples
	Float_t pedStdi; // pedestal standard deviation over all samples
	Int_t chan; // channel number
	Int_t fee; // fee number
	Int_t module; // module number (e.g. R1, R2, R3)
	Int_t slot; // slot number
	outputTree->Branch("isAlive",&isAlive,"isAlive/I"); 
	outputTree->Branch("pedMean",&pedMean,"pedMean/F");
     	outputTree->Branch("pedStdi",&pedStdi,"pedStdi/F");
        outputTree->Branch("chan",&chan,"chan/I");
        outputTree->Branch("fee",&fee,"fee/I");
        outputTree->Branch("module",&module,"module/I");
	outputTree->Branch("slot",&slot,"slot/I");

    

    string dChan = "Dead Channels (Sector "+sectorNumber+", Run "+runNumber+");Channel Number + 256*FEE Number;Instances"; 

    //TH1F *ave_adc = new TH1F("ave_adc","Channel Pedestals;Channel Number + 256*FEE Number;Pedestal ADC",256*26,-0.5,255.5+256*25);
    TH1F *dead_channel = new TH1F("dead_channel",dChan.c_str(),256*26,-0.5,255.5+256*25);
    //TH1F *rms_adc = new TH1F("rms_adc","Pedestal RMS;Channel Number + 256*FEE Number;Pedestal RMS",256*26,-0.5,255.5+256*25);
  
    std::vector<float> channels;
    std::vector<float> pedestal;
    std::vector<float> std_data;

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

    //Int_t channels[256];
    while(myReader.Next())
        {
            //cout<<"fee="<<*fee<<"Channel="<<*Channel<<endl;
            if(*nSamples_in==0) 
	    {
		dead_channel->Fill(*Channel_in+ 256*(*fee_in)); 
		alive_array_fee_channel[*fee_in][*Channel_in]=0; 
		continue;
	    }
            
            bool dead = false;
            for(int adc_sam_no=0;adc_sam_no<*nSamples_in;adc_sam_no++){
                //if (*fee_in == 25) channels[*Channel_in] += adcSamples_in[adc_sam_no];
                if (adcSamples_in[adc_sam_no] == 0 || TMath::IsNaN(float(adcSamples_in[adc_sam_no]))) {
			dead_channel->Fill(*Channel_in+ 256*(*fee_in));dead=true;break;
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

    for(int fee_no=0;fee_no<26;fee_no++){
          for(int channel_no=0;channel_no<256;channel_no++){
            if(counts_adc_fee_channel[fee_no][channel_no] != 0.0) 
            {
                float temp1 = ave_adc_fee_channel[fee_no][channel_no]/counts_adc_fee_channel[fee_no][channel_no];
                float temp2 = std_adc_fee_channel[fee_no][channel_no]/counts_adc_fee_channel[fee_no][channel_no];
                ave_adc_fee_channel[fee_no][channel_no] = temp1;
                std_adc_fee_channel[fee_no][channel_no] = temp2;
            }
            else
	    {
 		dead_channel->Fill(channel_no+ 256*(fee_no));
                //std::cout << "FEE: " << fee_no << ", channel: " << channel_no << std::endl;
                ave_adc_fee_channel[fee_no][channel_no] = 0.0;
                std_adc_fee_channel[fee_no][channel_no] = 0.0;
		alive_array_fee_channel[fee_no][channel_no]=0;
            }
            
	    if(sqrt(std_adc_fee_channel[fee_no][channel_no] - pow(ave_adc_fee_channel[fee_no][channel_no],2)) > 4)
  	    {
              outdata << "Sector " << q << ", FEE " << fee_no << ", Channel " << channel_no << std::endl;
 	    }
            if(ave_adc_fee_channel[fee_no][channel_no] > 200)
            {
 		dead_channel->Fill(channel_no+ 256*(fee_no));
  		std::cout << "(Large Pedestal pedestal>200) FEE: " << fee_no << ", channel: " << channel_no << ", " << ave_adc_fee_channel[fee_no][channel_no] << std::endl;
		alive_array_fee_channel[fee_no][channel_no]=0;
            }
            else if(ave_adc_fee_channel[fee_no][channel_no] < 10)
            {
 		dead_channel->Fill(channel_no+ 256*(fee_no));
  		std::cout << "(Small Pedestal pedestal<10) FEE: " << fee_no << ", channel: " << channel_no << ", " << ave_adc_fee_channel[fee_no][channel_no] << std::endl;
		alive_array_fee_channel[fee_no][channel_no]=0;
            }
            channels.push_back(channel_no+256.*(fee_no));
            pedestal.push_back(ave_adc_fee_channel[fee_no][channel_no]);
            std_data.push_back(sqrt(std_adc_fee_channel[fee_no][channel_no] - pow(ave_adc_fee_channel[fee_no][channel_no],2)));
            //ave_adc->SetBinContent(channel_no*(fee_no+1),ave_adc_fee_channel[fee_no][channel_no]);
            //std::cout << ave_adc_fee_channel[fee_no][channel_no] << std::endl;

	    pedMean=ave_adc_fee_channel[fee_no][channel_no];
	    pedStdi=sqrt(std_adc_fee_channel[fee_no][channel_no] - pow(ave_adc_fee_channel[fee_no][channel_no],2));
	    isAlive=alive_array_fee_channel[fee_no][channel_no];
            //if(*nSamples_in==0 || pedMean_o>200 || pedMean_o<10){isAlive=false};
	    fee=fee_no;
	    module=mod_arr[fee_no];
	    slot=slot_arr[fee_no];
	    outputTree->Fill();
	    

          }
      }
      //for (int l = 0; l < 256; l++)
      //{
      //    std::cout << "Channel " << l << ": " << ave_adc_fee_channel[25][l] << std::endl;
      //}

      TGraph *ave_adc = new TGraph(256*26,channels.data(),pedestal.data());
      TGraph *rms_adc = new TGraph(256*26,channels.data(),std_data.data());
      
      TCanvas *c1 = new TCanvas("c1","c1",2000,1200);
      c1->Divide(2,2);
      c1->cd(1);
      auto tl = new TLine(); tl->SetLineColor(kRed);
      dead_channel->Draw();
      for (int k = 0; k < 26; k++)
      {
        tl->DrawLine(256*(k+1),0,256*(k+1),dead_channel->GetMaximum());
        TLatex *pt = new TLatex(25+k*256,dead_channel->GetMaximum()*5/6,std::to_string(k).c_str());
        pt->SetTextSize(0.04);
   	pt->Draw();
        //TLatex *pt2 = new TLatex(25+k*256,dead_channel->GetMaximum()*7.5/10,mod[k].c_str());
        //pt2->SetTextSize(0.03);
   	//pt2->Draw();
        //TLatex *pt3 = new TLatex(25+k*256,dead_channel->GetMaximum()*4/6,slot[k].c_str());
        //pt3->SetTextSize(0.04);
   	//pt3->Draw();
      }
      c1->cd(2);
      ave_adc->SetTitle("Mean Pedestal Value;Channel Number + 256*FEE Number;Mean Pedestal ADC");
      ave_adc->SetMarkerStyle(20);
      ave_adc->SetMarkerColor(kBlue);
      ave_adc->SetMarkerSize(0.25);
      ave_adc->GetXaxis()->SetRangeUser(0,256*26);
      ave_adc->GetYaxis()->SetRangeUser(0,150);
      ave_adc->Draw("AP");
      for (int k = 0; k < 26; k++)
      {
        tl->DrawLine(256*(k+1),0,256*(k+1),150);
        TLatex *pt = new TLatex(25+k*256,135,std::to_string(k).c_str());
        pt->SetTextSize(0.04);
   	pt->Draw();
	/*
        //TLatex *pt2 = new TLatex(25+k*256,125,mod[k].c_str());
        //pt2->SetTextSize(0.03);
   	//pt2->Draw();
        TLatex *pt3 = new TLatex(25+k*256,115,slot[k].c_str());
        pt3->SetTextSize(0.04);
   	pt3->Draw();
      */
      }
      c1->cd(3);
      rms_adc->SetTitle("Pedestal Standard Deviation;Channel Number + 256*FEE Number;Pedestal Standard Deviation");
      rms_adc->SetMarkerStyle(20);
      rms_adc->SetMarkerColor(kBlue);
      rms_adc->SetMarkerSize(0.25);
      rms_adc->GetXaxis()->SetRangeUser(0,256*26);
      rms_adc->GetYaxis()->SetRangeUser(0,15);
      rms_adc->Draw("AP");
      for (int k = 0; k < 26; k++)
      {
        tl->DrawLine(256*(k+1),0,256*(k+1),15);
        TLatex *pt = new TLatex(25+k*256,13.5,std::to_string(k).c_str());
        pt->SetTextSize(0.04);
   	pt->Draw();
        //TLatex *pt2 = new TLatex(25+k*256,12.5,mod[k].c_str());
        //pt2->SetTextSize(0.03);
   	//pt2->Draw();
        //TLatex *pt3 = new TLatex(25+k*256,11.5,slot[k].c_str());
        //pt3->SetTextSize(0.04);
   	//pt3->Draw();
      }
 
      string saveName = "images/S"+sectorNumber+"-"+runNumber+".png";
      c1->SaveAs(saveName.c_str());
      c1->Close();
  outputfile->cd();
  outputTree->Write();
  outputfile->Close();
  }
  outdata.close();
}
