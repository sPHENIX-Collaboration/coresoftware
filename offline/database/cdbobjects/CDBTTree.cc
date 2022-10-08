#include "CDBTTree.h"

#include <phool/phool.h>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <climits>
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
  std::string fieldname = "F" + name;
  if (m_Locked[MultipleEntries])
  {
    std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
    std::cout << "That does not work, restructure your code" << std::endl;
    gSystem->Exit(1);
  }
  m_FloatEntryMap[channel].insert(std::make_pair(fieldname,value));
  m_EntryNameSet.insert(name);
}

void CDBTTree::SetDoubleValue(int channel, const std::string &name, double value)
{
  if (name == "ID")
  {
    std::cout << "Sorry ID is reserved as fieldname, pick anything else" << std::endl;
    gSystem->Exit(1);
  }
  std::string fieldname = "D" + name;
  if (m_Locked[MultipleEntries])
  {
    std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
    std::cout << "That does not work, restructure your code" << std::endl;
    gSystem->Exit(1);
  }
  m_DoubleEntryMap[channel].insert(std::make_pair(fieldname,value));
  m_EntryNameSet.insert(name);
}

void CDBTTree::SetIntValue(int channel, const std::string &name, int value)
{
  if (name == "ID")
  {
    std::cout << "Sorry ID is reserved as fieldname, pick anything else" << std::endl;
    gSystem->Exit(1);
  }
  std::string fieldname = "I" + name;
  if (m_Locked[MultipleEntries])
  {
    std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
    std::cout << "That does not work, restructure your code" << std::endl;
    gSystem->Exit(1);
  }
  m_IntEntryMap[channel].insert(std::make_pair(fieldname,value));
  m_EntryNameSet.insert(name);
}

void CDBTTree::SetUInt64Value(int channel, const std::string &name, uint64_t value)
{
  if (name == "ID")
  {
    std::cout << "Sorry ID is reserved as fieldname, pick anything else" << std::endl;
    gSystem->Exit(1);
  }
  std::string fieldname = "g" + name;
  if (m_Locked[MultipleEntries])
  {
    std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
    std::cout << "That does not work, restructure your code" << std::endl;
    gSystem->Exit(1);
  }
  m_UInt64EntryMap[channel].insert(std::make_pair(fieldname,value));
  m_EntryNameSet.insert(name);
}

void CDBTTree::Commit()
{
  m_Locked[MultipleEntries] = true;
// create ntuple string
  m_TTree[MultipleEntries] = new TTree(m_TTreeName[MultipleEntries].c_str(),m_TTreeName[MultipleEntries].c_str());
  std::set<int> id_set;
  std::map<std::string, float> floatmap;
  std::map<std::string, double> doublemap;
  std::map<std::string, int> intmap;
  std::map<std::string, uint64_t> uint64map;
  std::map<std::string, unsigned int> index;
  intmap.insert(std::make_pair("IID",INT_MIN));
  for (auto &f_entry: m_FloatEntryMap)
  {
    id_set.insert(f_entry.first);
    for (auto &f_val: f_entry.second)
    {
      floatmap.insert(std::make_pair(f_val.first,NAN));
    }
  }
  for (auto &f_val : floatmap)
  {
    std::string fielddescriptor =  f_val.first + "/F";
    m_TTree[MultipleEntries]->Branch(f_val.first.c_str(),&f_val.second,fielddescriptor.c_str());
  }

  for (auto &f_entry: m_DoubleEntryMap)
  {
    id_set.insert(f_entry.first);
    for (auto &f_val: f_entry.second)
    {
      doublemap.insert(std::make_pair(f_val.first,NAN));
    }
  }
  for (auto &f_val : doublemap)
  {
    std::string fielddescriptor =  f_val.first + "/D";
    m_TTree[MultipleEntries]->Branch(f_val.first.c_str(),&f_val.second,fielddescriptor.c_str());
  }

  for (auto &i_entry: m_IntEntryMap)
  {
    id_set.insert(i_entry.first);
    for (auto &i_val: i_entry.second)
    {
      intmap.insert(std::make_pair(i_val.first,INT_MIN));
    }
  }
  for (auto &i_val : intmap)
  {
    std::string fielddescriptor =  i_val.first + "/I";
    m_TTree[MultipleEntries]->Branch(i_val.first.c_str(),&i_val.second,fielddescriptor.c_str());
  }

  for (auto &i_entry: m_UInt64EntryMap)
  {
    id_set.insert(i_entry.first);
    for (auto &i_val: i_entry.second)
    {
      uint64map.insert(std::make_pair(i_val.first,UINT64_MAX));
    }
  }
  for (auto &i_val : uint64map)
  {
    std::string fielddescriptor =  i_val.first + "/g";
    m_TTree[MultipleEntries]->Branch(i_val.first.c_str(),&i_val.second,fielddescriptor.c_str());
  }
// fill ttree
  for (auto ids : id_set)
  {
    intmap["IID"] = ids;
    auto fmapiter = m_FloatEntryMap.find(ids);
    if (fmapiter !=  m_FloatEntryMap.end())
    {
      for (auto &f_val : fmapiter->second)
      {
	floatmap[f_val.first] = f_val.second;
      }
    }
    auto dmapiter = m_DoubleEntryMap.find(ids);
    if (dmapiter !=  m_DoubleEntryMap.end())
    {
      for (auto &d_val : dmapiter->second)
      {
	doublemap[d_val.first] = d_val.second;
      }
    }
    auto imapiter = m_IntEntryMap.find(ids);
    if (imapiter != m_IntEntryMap.end())
    {
      for (auto &i_val : imapiter->second)
      {
	intmap[i_val.first] = i_val.second;
      }
    }
    auto uint64mapiter = m_UInt64EntryMap.find(ids);
    if (uint64mapiter != m_UInt64EntryMap.end())
    {
      for (auto &uint64_val : uint64mapiter->second)
      {
	uint64map[uint64_val.first] = uint64_val.second;
      }
    }
    m_TTree[MultipleEntries]->Fill();
    for (auto &f_val: floatmap)
    {
      f_val.second = NAN;
    }
    for (auto &f_val: doublemap)
    {
      f_val.second = NAN;
    }
    for (auto &i_val: intmap)
    {
      i_val.second = INT_MIN;
    }
    for (auto &i_val: uint64map)
    {
      i_val.second = UINT64_MAX;
    }
  }
  return;
}

void CDBTTree::SetSingleFloatValue(const std::string &name, float value)
{
  std::string fieldname = "F" + name;
  if (m_SingleFloatEntryMap.find(fieldname) == m_SingleFloatEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleFloatEntryMap.insert(std::make_pair(fieldname, value));
    return;
  }
  m_SingleFloatEntryMap[fieldname] = value;
}

void CDBTTree::SetSingleDoubleValue(const std::string &name, double value)
{
  std::string fieldname = "D" + name;
  if (m_SingleDoubleEntryMap.find(fieldname) == m_SingleDoubleEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleDoubleEntryMap.insert(std::make_pair(fieldname, value));
    return;
  }
  m_SingleDoubleEntryMap[fieldname] = value;
}

void CDBTTree::SetSingleIntValue(const std::string &name, int value)
{
  std::string fieldname = "I" + name;
  if (m_SingleIntEntryMap.find(fieldname) == m_SingleIntEntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleIntEntryMap.insert(std::make_pair(fieldname, value));
    return;
  }
  m_SingleIntEntryMap[fieldname] = value;
}

void CDBTTree::SetSingleUInt64Value(const std::string &name, uint64_t value)
{
  std::string fieldname = "g" + name;
  if (m_SingleUInt64EntryMap.find(fieldname) == m_SingleUInt64EntryMap.end())
  {
    if (m_Locked[SingleEntries])
    {
      std::cout << "Trying to add field " << name << " after another entry was committed" << std::endl;
      std::cout << "That does not work, restructure your code" << std::endl;
      gSystem->Exit(1);
    }
    m_SingleUInt64EntryMap.insert(std::make_pair(fieldname, value));
    return;
  }
  m_SingleUInt64EntryMap[fieldname] = value;
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
      std::string tmpstring = calibs.first;
      tmpstring.erase(0);
	std::cout << "name " << tmpstring << " value: " << calibs.second << std::endl;
      }
    }
  }
  if (!m_SingleFloatEntryMap.empty())
  {
    std::cout << "Number of single float fields: " << m_SingleFloatEntryMap.size() << std::endl;
    for (auto &field : m_SingleFloatEntryMap)
    {
      std::string tmpstring = calibs.first;
      tmpstring.erase(0);
      std::cout << tmpstring << " value " << field.second << std::endl;
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
      std::string tmpstring = field.first;
      tmpstring.erase(0);
      std::cout << tmpstring << " value " << field.second << std::endl;
    }
    std::cout.copyfmt(oldState);
  }
  if (!m_SingleIntEntryMap.empty())
  {
    std::cout << "Number of single int fields: " << m_SingleIntEntryMap.size() << std::endl;
    for (auto &field : m_SingleIntEntryMap)
    {
      std::string tmpstring = field.first;
      tmpstring.erase(0);
      std::cout << tmpstring << " value " << field.second << std::endl;
    }
  }
  if (!m_SingleUInt64EntryMap.empty())
  {
    std::cout << "Number of single uint64 fields: " << m_SingleUInt64EntryMap.size() << std::endl;
    for (auto &field : m_SingleUInt64EntryMap)
    {
      std::string tmpstring = field.first;
      tmpstring.erase(0);
      std::cout << tmpstring << " value " << field.second << std::endl;
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
  std::string fieldname = "F" + name;
  auto singleiter = m_SingleFloatEntryMap.find(fieldname);
  if (singleiter == m_SingleFloatEntryMap.end())
  {
    std::cout << "Could not find " << name << " in single calibrations" << std::endl;
    std::cout << "Existing values:" << std::endl;
    for (auto &eiter : m_SingleFloatEntryMap)
    {
      std::string tmpstring = eiter.first;
      tmpstring.erase(0);
      std::cout << "name : " << tmpstring << ", value " << eiter.second
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
  std::string fieldname = "F" + name;
  auto calibiter = channelmapiter->second.find(fieldname);
  if (calibiter == channelmapiter->second.end())
  {
    std::cout << "Could not find " << name << " among calibrations" << std::endl;
    return NAN;
  }
  return calibiter->second;
}
