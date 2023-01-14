#include "CDBHistos.h"

#include <phool/phool.h>

#include <TClass.h>       // for TClass
#include <TCollection.h>  // for TIter
#include <TDirectory.h>   // for TDirectoryAtomicAdapter, TDirectory, gDirec...
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>    // for TList
#include <TObject.h>  // for TObject
#include <TROOT.h>
#include <TSystem.h>

#include <iostream>
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
  if (m_HistoMap.empty())
  {
    std::cout << PHWHERE << " no histograms to be saved " << std::endl;
    return;
  }
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
  if (fin == nullptr)
  {
    std::cout << PHWHERE << " Could not open " << m_Filename << std::endl;
    return;
  }
  TList *list = fin->GetListOfKeys();
  if (!list)
  {
    std::cout << PHWHERE << " No keys found in " << m_Filename << std::endl;
    fin->Close();
    gROOT->cd(currdir.c_str());  // restore previous directory
    return;
  }
  TIter next(list);
  TKey *key;
  TObject *obj;
  TH1 *h1 = nullptr;
  while ((key = (TKey *) next()))
  {
    obj = key->ReadObj();
    if ((obj->InheritsFrom("TH1")))
    {
      fin->GetObject(obj->GetName(), h1);
      h1->SetDirectory(nullptr);
      m_HistoMap.insert(std::make_pair(obj->GetName(), h1));
    }
  }
  fin->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBHistos::Print() const
{
  for (auto &iter : m_HistoMap)
  {
    std::cout << "histogram " << iter.first << ", type "
              << iter.second->IsA()->GetName() << std::endl;
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
  m_HistoMap.insert(std::make_pair(h1->GetName(), h1));
  return;
}

TH1 *CDBHistos::getHisto(const std::string &name)
{
  const auto iter = m_HistoMap.find(name);
  if (iter == m_HistoMap.end())
  {
    std::cout << PHWHERE << ": Histogram " << name << " not found" << std::endl;
    return nullptr;
  }
  return iter->second;
}
