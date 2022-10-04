#include "CDBTTree.h"

#include <phool/phool.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <iostream>
#include <limits>

CDBTTree::CDBTTree(const std::string &fname)
  : m_Filename(fname)
{}

CDBTTree::~CDBTTree()
{
  m_FloatEntryMap.clear();
  m_SingleFloatEntryMap.clear();
}

void CDBTTree::SetFloatValue(int channel, const std::string &name, float value)
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
  m_FloatEntryMap[channel].insert(std::make_pair(name,value));
  m_EntryNameSet.insert(name);
}

void CDBTTree::Commit()
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
  float *vals = new float[m_EntryNameSet.size()+1]; // +1 for the id
  std::fill(vals,vals+m_EntryNameSet.size()+1,NAN);
  if (m_TTree[MultipleEntries] == nullptr)
  {
//    m_TTree[MultipleEntries] = new TTree(m_TTreeName[MultipleEntries].c_str(),m_TTreeName[MultipleEntries].c_str(),fieldnames.c_str());
    m_TTree[MultipleEntries] = new TTree(m_TTreeName[MultipleEntries].c_str(),m_TTreeName[MultipleEntries].c_str());
  }
  for (auto &entry :  m_FloatEntryMap)
  {
    vals[(index.find("ID")->second)] = entry.first;
    for (auto &calib : entry.second)
    {
      vals[(index.find(calib.first)->second)] = calib.second;  
    }
//    m_TTree[MultipleEntries]->Fill(vals);
    std::fill(vals,vals+m_EntryNameSet.size()+1,NAN);
  }
  delete [] vals;
  return;
}

void CDBTTree::SetSingleFloatValue(const std::string &name, float value)
{
  if (m_SingleFloatEntryMap.find(name) == m_SingleFloatEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleFloatEntryMap.insert(std::make_pair(name, value));
    return;
  }
  m_SingleFloatEntryMap[name] = value;
}

void CDBTTree::SetSingleDoubleValue(const std::string &name, double value)
{
  if (m_SingleDoubleEntryMap.find(name) == m_SingleDoubleEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleDoubleEntryMap.insert(std::make_pair(name, value));
    return;
  }
  m_SingleDoubleEntryMap[name] = value;
}

void CDBTTree::SetSingleIntValue(const std::string &name, int value)
{
  if (m_SingleIntEntryMap.find(name) == m_SingleIntEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleIntEntryMap.insert(std::make_pair(name, value));
    return;
  }
  m_SingleIntEntryMap[name] = value;
}

void CDBTTree::SetSingleUInt64Value(const std::string &name, uint64_t value)
{
  if (m_SingleUInt64EntryMap.find(name) == m_SingleUInt64EntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleUInt64EntryMap.insert(std::make_pair(name, value));
    return;
  }
  m_SingleUInt64EntryMap[name] = value;
}

void CDBTTree::CommitSingle()
{
  m_Locked[SingleEntries] = true;
  m_TTree[SingleEntries] = new TTree(m_TTreeName[SingleEntries].c_str(),m_TTreeName[SingleEntries].c_str());
  for (auto &field :  m_SingleFloatEntryMap)
  {
    std::string fielddescriptor =  field.first + "/F";
    m_TTree[SingleEntries]->Branch(field.first.c_str(),&field.second,fielddescriptor.c_str());
  }
  for (auto &field :  m_SingleDoubleEntryMap)
  {
    std::string fielddescriptor =  field.first + "/D";
    m_TTree[SingleEntries]->Branch(field.first.c_str(),&field.second,fielddescriptor.c_str());
  }
  for (auto &field :  m_SingleIntEntryMap)
  {
    std::string fielddescriptor =  field.first + "/I";
    m_TTree[SingleEntries]->Branch(field.first.c_str(),&field.second,fielddescriptor.c_str());
  }
  for (auto &field :  m_SingleUInt64EntryMap)
  {
    std::string fielddescriptor =  field.first + "/g";
    m_TTree[SingleEntries]->Branch(field.first.c_str(),&field.second,fielddescriptor.c_str());
  }

  m_TTree[SingleEntries]->Fill();
  return;
}

void CDBTTree::Print()
{
  if (! m_FloatEntryMap.empty())
  {
    std::cout << "Number of float entries: " << m_FloatEntryMap.size() << std::endl;
    for (auto &field : m_FloatEntryMap)
    {
      std::cout << "ID: " << field.first << std::endl;
      for (auto &calibs : field.second)
      {
	std::cout << "name " << calibs.first << " value: " << calibs.second << std::endl;
      }
    }
  }
  if (!m_SingleFloatEntryMap.empty())
  {
    std::cout << "Number of single float fields: " << m_SingleFloatEntryMap.size() << std::endl;
    for (auto &field : m_SingleFloatEntryMap)
    {
      std::cout << field.first << " value " << field.second << std::endl;
    }
  }
  if (!m_SingleDoubleEntryMap.empty())
  {
    std::cout << "Number of single double fields: " << m_SingleDoubleEntryMap.size() << std::endl;
// some acrobatics to restore the old state of cout after changing the precision for double printout
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);
    std::cout.precision(std::numeric_limits< double >::max_digits10);
    for (auto &field : m_SingleDoubleEntryMap)
    {
      std::cout << field.first << " value " << field.second << std::endl;
    }
    std::cout.copyfmt(oldState);
  }
  if (!m_SingleIntEntryMap.empty())
  {
    std::cout << "Number of single int fields: " << m_SingleIntEntryMap.size() << std::endl;
    for (auto &field : m_SingleIntEntryMap)
    {
      std::cout << field.first << " value " << field.second << std::endl;
    }
  }
  if (!m_SingleUInt64EntryMap.empty())
  {
    std::cout << "Number of single uint64 fields: " << m_SingleUInt64EntryMap.size() << std::endl;
    for (auto &field : m_SingleUInt64EntryMap)
    {
      std::cout << field.first << " value " << field.second << std::endl;
    }
  }
}

void CDBTTree::WriteCDBTTree()
{
  TFile *f = TFile::Open(m_Filename.c_str(),"RECREATE");
  for (auto ntup : m_TTree)
  {
    if (ntup != nullptr)
    {
      ntup->Write();
    }
  }
  f->Close();
}

void CDBTTree::LoadCalibrations()
{
  TFile *f = TFile::Open(m_Filename.c_str());
  f->GetObject(m_TTreeName[SingleEntries].c_str(), m_TTree[SingleEntries]);
  f->GetObject(m_TTreeName[MultipleEntries].c_str(), m_TTree[MultipleEntries]);
  if (m_TTree[SingleEntries] != nullptr)
  {
    TObjArray *branches =  m_TTree[SingleEntries]->GetListOfBranches();
    TIter iter(branches);
    while(TBranch *thisbranch = static_cast<TBranch *> (iter.Next()))
    {
      // this convoluted expression returns the data type of a split branch
      std::string DataType = thisbranch->GetLeaf(thisbranch->GetName())->GetTypeName();
      if (DataType == "Float_t")
      {
	auto itermap =  m_SingleFloatEntryMap.insert(std::make_pair(thisbranch->GetName(),NAN));
	m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&(itermap.first)->second);
      }
      else if (DataType == "Double_t")
      {
	auto itermap =  m_SingleDoubleEntryMap.insert(std::make_pair(thisbranch->GetName(),NAN));
	m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&(itermap.first)->second);
      }
      else if (DataType == "Int_t")
      {
	auto itermap =  m_SingleIntEntryMap.insert(std::make_pair(thisbranch->GetName(),-99999));
	m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&(itermap.first)->second);
      }
      else if (DataType == "ULong_t")
      {
	auto itermap =  m_SingleUInt64EntryMap.insert(std::make_pair(thisbranch->GetName(),~0x0));
	m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&(itermap.first)->second);
      }
      else
      {
	std::cout << PHWHERE << " data type " << DataType
		  << " in " << m_TTree[SingleEntries]->GetName()
		  << " from " << f->GetName()
                  << " not implemented" << std::endl;
	gSystem->Exit(1);
      }
//      thisbranch->Print();
      // m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(),&val);
      // thisbranch->GetEntry(0);
      // m_SingleFloatEntryMap.insert(std::make_pair(thisbranch->GetName(), val));
    }
  }
m_TTree[SingleEntries]->GetEntry(0);
  if (m_TTree[MultipleEntries] != nullptr)
  {
    TObjArray *branches =  m_TTree[MultipleEntries]->GetListOfBranches();
    TIter iter(branches);
    std::map<std::string, float> branchmap;
    while(TBranch *thisbranch = static_cast<TBranch *> (iter.Next()))
    {
      float val;
      std::pair<std::string, float> valpair = std::make_pair(thisbranch->GetName(),val);
      branchmap.insert(valpair);
    }
    for (auto &biter : branchmap)
    {
      m_TTree[MultipleEntries]->SetBranchAddress(biter.first.c_str(), &biter.second);
    }
    for (int i=0; i<m_TTree[MultipleEntries]->GetEntries(); i++)
    {
      m_TTree[MultipleEntries]->GetEntry(i);
      int id = branchmap.find("ID")->second;
      m_FloatEntryMap.insert(std::make_pair(id, branchmap));
    }
  }
  f->Close();
}

float CDBTTree::GetSingleFloatValue(const std::string &name)
{  
  if ( m_SingleFloatEntryMap.empty())
  {
    LoadCalibrations();
  }
  auto singleiter = m_SingleFloatEntryMap.find(name);
  if (singleiter == m_SingleFloatEntryMap.end())
  {
    std::cout << "Could not find " << name << " in single calibrations" << std::endl;
    std::cout << "Existing values:" << std::endl;
    for (auto &eiter : m_SingleFloatEntryMap)
    {
      std::cout << "name : " << eiter.first << ", value " << eiter.second
		<< std::endl;
    }
  }
  return singleiter->second ;
}

float CDBTTree::GetFloatValue(int channel, const std::string &name)
{
  if ( m_FloatEntryMap.empty())
  {
    LoadCalibrations();
  }
  auto channelmapiter = m_FloatEntryMap.find(channel);
  if (channelmapiter == m_FloatEntryMap.end())
  {
    std::cout << "Could not find channel " << channel << " in calibrations" << std::endl;
    return NAN;
  }
  auto calibiter = channelmapiter->second.find(name);
  if (channelmapiter->second.find(name) == channelmapiter->second.end())
  {
    std::cout << "Could not find " << name << " among calibrations" << std::endl;
    return NAN;
  }
  return calibiter->second;
}
