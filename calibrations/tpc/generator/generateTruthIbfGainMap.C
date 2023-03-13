void generateTruthIbfGainMap(const char* adcFile, const char *adcName, const char* ionFile, const char* ibfName, const char* primName, const char* outputFile){
  
 //load the adc-per-bin data from the specified file.
  TFile* adcInputFile = TFile::Open(adcFile, "READ");
  TH3* hAdc=nullptr;
  adcInputFile->GetObject(adcName, hAdc);
  if (hAdc == nullptr)  {
    printf("ADC hist %s or file %s not found!\n",adcName, adcFile);
    return false;
  }
  
  //load the ions-per-bin data from the specified file
  TFile* ibfInputFile = TFile::Open(ionFile, "READ");
  TH2* hIbf=nullptr;
  TH2* hPrimaries=nullptr;
  ibfInputFile->GetObject(ibfName,hIbf);
  ibfInputFile->GetObject(primName,hPrimaries);
  if (hIbf == nullptr) {
    printf("IBF hist %s or file %s not found!\n",ibfName, ionFile);
    return false;
  }
  if (hPrimaries == nullptr) {
    printf("IBF hist %s IN file %s not found!\n",primName, ionFile);
    return false;
  }

  TFile *output=TFile::Open(outputFile,"RECREATE");
  //generate a histogram with the IBF + Primaries:
  TH3* hTotalCharge=nullptr;
  hTotalCharge=static_cast<TH3*>(hPrimaries->Clone("hTotalCharge"));
  TH2* hFlatTotal=(TH2*)hTotalCharge->Project3D("xy");
  TH2* hFlatIbf=(TH2*)hIbf->Project3D("xy");
  
  
  
  //sanity check that the bounds are the same

  TAxis* ax[3] = {nullptr, nullptr, nullptr};
  ax[0] = hPrimaries->GetXaxis();
  ax[1] = hPrimaries->GetYaxis();
  ax[2] = hPrimaries->GetZaxis();

  TAxis* ax2[3] = {nullptr, nullptr, nullptr};
  ax2[0] = hAdc->GetXaxis();
  ax2[1] = hAdc->GetYaxis();
  ax2[2] = hAdc->GetYaxis();


  int nbins[3];
  for (int i = 0; i < 3; i++)
  {
    nbins[i] = ax[i]->GetNbins();  //number of bins, not counting under and overflow.
    if (nbins[i]!=ax2[i]->GetNbins()){
      printf("Primaries and Adc bins are different in axis %d.  Failing\n",i);
      return;
    }
  }

  
  //project into 2D per rphi, and divide Ions by Adc to get the Ion/Adc gain.
  TH2* hIonGain[2],*hFlatAdc;
  hIonGain[0]=static_cast<TH2*>(hFlatTotal->Clone("hIonGain"));//with primaries
  hIonGain[1]=static_cast<TH2*>(hFlatIbf->Clone("hIbfGain"));//with only ibf
  hFlatAdc=static_cast<TH2*>(hAdc->Project3D("xy"));

  
  hIonGain[0]->Divide(hFlatAdc);
  hIonGain[1]->Divide(hFlatAdc);

  hIonGain[0]->Write();
  hIonGain[1]->Write();
  
  //
  return;
}
