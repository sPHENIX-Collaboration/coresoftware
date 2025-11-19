#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>

#include <format>
#include <iostream>

bool generateTruthIbfGainMap(const char* adcFile, const char *adcName, const char* ionFile, const char* ibfName, const char* primName, const char*  /*outputFile*/){
  
 //load the adc-per-bin data from the specified file.
  TFile* adcInputFile = TFile::Open(adcFile, "READ");
  TH3* hAdc=nullptr;
  adcInputFile->GetObject(adcName, hAdc);
  if (hAdc == nullptr)  {
    std::cout << "ADC hist " <<  adcName << " or file "
	      <<  adcFile << " not found!" << std::endl;
    return false;
  }
  
  //load the ions-per-bin data from the specified file
  TFile* ibfInputFile = TFile::Open(ionFile, "READ");
  TH3* hIbf=nullptr;
  TH3* hPrimaries=nullptr;
  ibfInputFile->GetObject(ibfName,hIbf);
  ibfInputFile->GetObject(primName,hPrimaries);
  if (hIbf == nullptr) {
    std::cout << "IBF hist " << ibfName << " or file "
	      << ionFile << " not found!" << std::endl;
    return false;
  }
  if (hPrimaries == nullptr) {
    std::cout << "IBF hist " <<  primName << " IN file "
	      << ionFile << " not found!" << std::endl;
    return false;
  }

//  TFile *output=TFile::Open(outputFile,"RECREATE");
  //generate 2D histograms with the IBF, and IBF + Primaries, per side:

  const int nSections=3;
  TH2 *hFlatTotal[nSections];
  TH2 *hFlatIbf[nSections]; //flat charge sums in north and south
  int zbins;
  int neg;
  int cm;
  int pos;
  zbins=hAdc->GetZaxis()->GetNbins();
  neg=1;
  cm=zbins/2;
  pos=zbins;

  TString suffix[]={"","_negz","_posz"};
  int low[]={neg,neg,cm+1};
  int high[]={pos,cm,pos};

  for (int i=0;i<nSections;i++){
    hPrimaries->GetZaxis()->SetRange(low[i],high[i]);
    hFlatTotal[i]=(TH2*)hPrimaries->Project3D(std::format("xy{}",suffix[i].Data()).c_str());
    hIbf->GetZaxis()->SetRange(low[i],high[i]);
    hFlatIbf[i]=(TH2*)hIbf->Project3D(std::format("xy{}",suffix[i].Data()).c_str());
    hFlatTotal[i]->Add(hFlatIbf[i]);
  }
  
  
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
  for (int i = 0; i < 3; i++)//note these are dimensions, not sections.
  {
    nbins[i] = ax[i]->GetNbins();  //number of bins, not counting under and overflow.
    if (nbins[i]!=ax2[i]->GetNbins()){
      std::cout << "Primaries and Adc bins are different in axis " << i
		<< ". (" << nbins[i] << " vs "<< ax2[i]->GetNbins() << ").  Aborting unless in z." << std::endl;
      if (i!=2)
      {
	return false;
      }
      std::cout << "this is a z difference, which we can recover." << std::endl;
    }
  }

  
  //project into 2D per rphi, and divide Ions by Adc to get the Ion/Adc gain.
  TH2 * hIonGain[nSections];
  TH2 *hIbfGain[nSections];
  TH2 *hFlatAdc[nSections];

  for (int i=0;i<nSections;i++){
    hIonGain[i]=static_cast<TH2*>(hFlatTotal[i]->Clone(std::format("hIonGain{}",suffix[i].Data()).c_str()));//with primaries
    hIbfGain[i]=static_cast<TH2*>(hFlatIbf[i]->Clone(std::format("hIbfGain{}",suffix[i].Data()).c_str()));//with only ibf
    hAdc->GetZaxis()->SetRange(low[i],high[i]);
    hFlatAdc[i]=static_cast<TH2*>(hAdc->Project3D(std::format("xy{}",suffix[i].Data()).c_str()));

    hIonGain[i]->Divide(hFlatAdc[i]);
    hIbfGain[i]->Divide(hFlatAdc[i]);

    hIonGain[i]->Write();
    hIbfGain[i]->Write();
  }
  
  //
  return true;
}
