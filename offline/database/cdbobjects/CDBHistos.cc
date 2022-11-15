#include "CDBHistos.h"

#include <phool/phool.h>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

#include <climits>
#include <cmath>    // for NAN, isfinite
#include <cstdint>  // for uint64_t, UINT64_MAX
#include <iostream>
#include <set>      // for set
#include <utility>  // for pair, make_pair

CDBHistos::CDBHistos(const std::string &fname)
  : m_Filename(fname)
{
}

CDBHistos::~CDBHistos()
{
  for (auto &iter : m_HistoMap)
  {
    delete iter.second;
  }
  m_HistoMap.clear();
}



void CDBHistos::WriteCDBHistos()
{
  std::string currdir = gDirectory->GetPath();
  TFile *f = TFile::Open(m_Filename.c_str(), "RECREATE");
  for (auto &iter : m_HistoMap)
  {
    iter.second->Write();
  }
  f->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBHistos::LoadCalibrations()
{
  std::string currdir = gDirectory->GetPath();

  TFile *fin = TFile::Open(m_Filename.c_str());
TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
  std::cout << "got keys" << std::endl;
  TH1 *h1 = nullptr;
  while ( (key = (TKey*)next()) ) 
  {
    obj = key->ReadObj() ;
    std::cout << "found key " << obj->GetName() << std::endl;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
	    || (obj->InheritsFrom("TH1"))
      ) 
    {
      fin->GetObject(obj->GetName(),h1);
h1->SetDirectory(nullptr);
      printf("<W> Object %s is 1D or 2D histogram : "
             ,obj->GetName()) ;
      m_HistoMap.insert(std::make_pair(obj->GetName(),h1));
    }
//    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());
  }
  fin->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBHistos::Print() const
{
  for (auto &iter : m_HistoMap)
  {
    std::cout << "histogram " << iter.first << std::endl;
    iter.second->Draw();
    sleep(10);
  }
  return;
}

void CDBHistos::registerHisto(TH1 *h1)
{
  const auto iter = m_HistoMap.find(h1->GetName());
  if (iter != m_HistoMap.end())
  {
    std::cout << PHWHERE << " Histogram " << h1->GetName() << " already registered, use a different name and try again" << std::endl;
    gSystem->Exit(1);
  }
  m_HistoMap.insert(std::make_pair(h1->GetName(),h1));
  return;
}
