#include "CDBNtuple.h"

#include <TFile.h>
#include <TNtuple.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <iostream>

CDBNtuple::CDBNtuple(const std::string &fname)
  : m_Filename(fname)
{
}

void CDBNtuple::SetValue(int channel, const std::string &name, float value)
{
  if (name == "ID")
  {
    std::cout << "Sorry ID is reserved as fieldname, pick anything else" << std::endl;
    gSystem->Exit(1);
  }
    if (m_Locked[MultipleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_EntryMap[channel].insert(std::make_pair(name,value));
    m_EntryNameSet.insert(name);
}

void CDBNtuple::Commit()
{
  m_Locked[MultipleEntries] = true;
// create ntuple string
  std::map<std::string, unsigned int> index;
  index.insert(std::make_pair("ID",0));
  std::string fieldnames = "ID";
  unsigned int i = 0;
  for (auto &field :  m_EntryNameSet)
  {
    fieldnames += ":";
    fieldnames += field;
    i++;
    index.insert(std::make_pair(field,i));
  }
  std::cout << "fieldstring: " << fieldnames << std::endl;
  float *vals = new float[m_EntryNameSet.size()+1]; // +1 for the id
  std::fill(vals,vals+m_EntryNameSet.size()+1,NAN);
  if (m_TNtuple[MultipleEntries] == nullptr)
  {
m_TNtuple[MultipleEntries] = new TNtuple(m_TNtupleName[MultipleEntries].c_str(),m_TNtupleName[MultipleEntries].c_str(),fieldnames.c_str());
  }
  for (auto &entry :  m_EntryMap)
  {
    std::cout << "entrymap entry " <<  entry.first 
	      << " index " << index.find("ID")->second << std::endl;
    vals[(index.find("ID")->second)] = entry.first;
    for (auto &calib : entry.second)
    {
      std::cout << "Committing " << calib.first << " at index " << index.find(calib.first)->second
		<< " value " <<  calib.second << std::endl; 
      vals[(index.find(calib.first)->second)] = calib.second;  
	 }
m_TNtuple[MultipleEntries]->Fill(vals);
  // for (size_t i = 0; i < m_EntryNameSet.size()+1; i++)
  // {
  //   vals[i] = NAN;
  // }
  std::fill(vals,vals+m_EntryNameSet.size()+1,NAN);
  }
  delete [] vals;
//  std::cout << calibid << std::endl;

  return;
}

void CDBNtuple::SetSingleValue(const std::string &name, float value)
{
  if (m_SingleEntryMap.find(name) == m_SingleEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleEntryMap.insert(std::make_pair(name, value));
    return;
  }
  m_SingleEntryMap[name] = value;
  std::cout << value << std::endl;
}

void CDBNtuple::CommitSingle()
{
  m_Locked[SingleEntries] = true;
  float *vals = new float[m_SingleEntryMap.size()];
  std::string fieldnames;
  size_t i=0;
  for (auto &field :  m_SingleEntryMap)
  {
    fieldnames += field.first;
    vals[i] =   field.second;  
    i++;
    if (i < m_SingleEntryMap.size())
    {
      fieldnames += ':';
    }
    field.second = NAN;
  }
  std::cout << "fields: " << fieldnames << std::endl;
  m_TNtuple[SingleEntries] = new TNtuple(m_TNtupleName[SingleEntries].c_str(),m_TNtupleName[SingleEntries].c_str(),fieldnames.c_str());
  m_TNtuple[SingleEntries]->Fill(vals);
  delete [] vals;
  std::cout << "commit single" << std::endl;
  return;
}

void CDBNtuple::Print()
{
  std::cout << "Number of fields: " << m_EntryMap.size() << std::endl;
  for (auto &field : m_EntryMap)
  {
    std::cout << field.first << std::endl;
    for (auto &calibs : field.second)
    {
      std::cout << "name " << calibs.first << " value: " << calibs.second << std::endl;
    }
  }
  std::cout << "Number of single fields: " << m_SingleEntryMap.size() << std::endl;
  for (auto &field : m_SingleEntryMap)
  {
    std::cout << field.first << " value " << field.second << std::endl;
  }
}

void CDBNtuple::WriteCDBNtuple()
{
  TFile *f = TFile::Open(m_Filename.c_str(),"RECREATE");
  for (auto ntup : m_TNtuple)
  {
    if (ntup != nullptr)
    {
      ntup->Write();
    }
  }
  f->Close();
}

void CDBNtuple::LoadCalibrations()
{
  TFile *f = TFile::Open(m_Filename.c_str());
  f->GetObject(m_TNtupleName[SingleEntries].c_str(), m_TNtuple[SingleEntries]);
  f->GetObject(m_TNtupleName[MultipleEntries].c_str(), m_TNtuple[MultipleEntries]);
  if (m_TNtuple[SingleEntries] != nullptr)
  {
    TObjArray *branches =  m_TNtuple[SingleEntries]->GetListOfBranches();
    TIter iter(branches); //= branches->MakeIterator();
    while(TBranch *thisbranch = static_cast<TBranch *> (iter.Next()))
    {
      float val;
      m_TNtuple[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&val);
      thisbranch->GetEntry(0);
      m_SingleEntryMap.insert(std::make_pair(thisbranch->GetName(), val));
    }
  }
  f->Close();
}

float CDBNtuple::GetSingleValue(const std::string &name)
{  
  if ( m_SingleEntryMap.empty())
 {
  LoadCalibrations();
}
  if (m_SingleEntryMap.find(name) == m_SingleEntryMap.end())
  {
    std::cout << "Could not find " << name << " in single calibrations" << std::endl;
    std::cout << "Existing values:" << std::endl;
  for (auto &eiter : m_SingleEntryMap)
  {
    std::cout << "name : " << eiter.first << ", value " << eiter.second 
	      << std::endl;
  }
  }
  return NAN;
}

