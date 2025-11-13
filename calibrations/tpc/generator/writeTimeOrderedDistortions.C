//basic framework to read and compare events in a time ordered distortion file
//these files contain a TTree with xingnum, and three distortion map branches

#include <TCanvas.h>
#include <TFile.h>
#include <TFileCollection.h>
#include <TFileInfo.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THashList.h>
#include <TTree.h>

#include <format>
#include <iostream>

void writeTimeOrderedDistortions(bool subtractFirst=false, const std::string &filename="/sphenix/user/rcorliss/distortion_maps/2022.07/TimeOrderedDistortions.root", const std::string &inputpattern="/sphenix/user/rcorliss/distortion_maps/2022.07/*.distortion_map.hist.root", bool debug=false){

  TFile *treefile=TFile::Open(filename.c_str(),"RECREATE");
  //TTree *tree=(TTree*)(file->Get("TimeDists"));

  const int nhists=8;
  TH3F *basehist[nhists];
  TH3F *temphist[nhists];
  int xingnum=0;
  //note that we match the histogram RPhi to the branch P, and the histogram P to the branch Phi, to maintain backwards compatibility with older sets... for now.
  std::string branchname[]={"hIntDistortionP_negz",
			    "hIntDistortionPhi_negz",
			    "hIntDistortionR_negz",
			    "hIntDistortionZ_negz",
			    "hIntDistortionP_posz",
			    "hIntDistortionPhi_posz",
			    "hIntDistortionR_posz",
			    "hIntDistortionZ_posz"};
  std::string histname[]={"hIntDistortionRPhi_negz",
			  "hIntDistortionP_negz",
			  "hIntDistortionR_negz",
			  "hIntDistortionZ_negz",
			  "hIntDistortionRPhi_posz",
			  "hIntDistortionP_posz",
			  "hIntDistortionR_posz",
			  "hIntDistortionZ_posz"};
  TTree *tree=new TTree("TimeDists", "TimeDists");
  tree->Branch("xingnum",&xingnum);
  for (int i=0;i<nhists;i++){
    temphist[i]=new TH3F(std::format("temphist{}",i).c_str(),std::format("temphist{}",i).c_str(),10,0,10,20,0,20,30,0,30);
    basehist[i]=new TH3F(std::format("basehist{}",i).c_str(),std::format("basehist{}",i).c_str(),10,0,10,20,0,20,30,0,30);
    tree->Branch(branchname[i].c_str(),&(temphist[i]));
  }
  if (debug) {  std::cout << "histograms built and branched." << std::endl;
}

  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern.c_str());
  filelist->Print();
  std::cout << "found " << filelist->GetNFiles() << " files like: " << ((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile() << std::endl;

  //return;
  TFile *infile;
  bool fileIsValid=true;
  bool isFirst=true;
  int nMaps=0;
  

  for (int i=0;i<filelist->GetNFiles() && i<100;i++){
    //for each file, find all histograms in that file.
    infile=TFile::Open(((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl(),"READ");//gross.
    fileIsValid=true;
    if (debug)
    {
      std::cout << "=====> Trying File " << i << std::endl;
    }
    if (!infile->IsOpen())
    {
      std::cout << "=====> File " << i << " is NOT openable <=======" << std::endl;
      continue; //file didn't open right.  move on to the next one.
    }
    //TList *keys=infile->GetListOfKeys();
    for (int j=0;j<nhists;j++){
      temphist[j]=nullptr;
      temphist[j]=infile->Get<TH3F>(histname[j].c_str());
      if (!temphist[j]){
	fileIsValid=false; //histogram doesn't exist.  don't bother loading the other hists.
	break;
      }
      int nbins=temphist[j]->GetNcells();
        if (debug)
	{
	  std::cout << "=======> \"" << histname[j] << "\" has " << nbins << " cells" << std::endl;
	}

    }
    if (!fileIsValid)
    {
      infile->Close();
      std::cout << "=====> File " << i << " is NOT valid <=======" << std::endl;

      continue; //didn't get all our hists.  move on to the next file.
    }
    if (debug)
    {
      std::cout << "=====> File " << i << " is valid" << std::endl;
    }
    xingnum=i;//temporary fix to paste something in there.
    if(subtractFirst){
      if (isFirst){
	for (int j=0;j<nhists;j++){
	  treefile->cd();
	  basehist[j]=(TH3F*)(temphist[j]->Clone());
	  //temphist[j]->Copy(*(basehist[j]));
	}
	isFirst=false;
      }
      for (int j=0;j<nhists;j++){
	std::cout << "ptr j=" << j << ":  b:" << basehist[j] << "\tt:" << temphist[j] << std::endl;
	int nbins=temphist[j]->GetNcells();
	for (int k=0;k<nbins;k++){
	  double b=basehist[j]->GetBinContent(k);
	  double t=temphist[j]->GetBinContent(k);
	  double diff=t-b;
	  temphist[j]->SetBinContent(k,diff);
	  //temphist[j]->Add(basehist[j],-1);
	}
	if (debug)
	{
	  std::cout << "=======> \"" << histname[j] << "\" has " << nbins << " cells when writing diff" << std::endl;
	}

      }
    }
    
    tree->Fill();
    nMaps++;
    infile->Close();
  }
  std::cout << "finished tree " << filename << " has nMaps=" << nMaps << std::endl;
      
  treefile->cd();
  tree->Write();
  treefile->Close();
  return;
}
