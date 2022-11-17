#include "CDBTTree.h"

#include <phool/phool.h>

#include <TBranch.h>      // for TBranch
#include <TCollection.h>  // for TIter
#include <TDirectory.h>   // for TDirectoryAtomicAdapter, TDirectory, gDirec...
#include <TFile.h>
#include <TLeaf.h>      // for TLeaf
#include <TObjArray.h>  // for TObjArray
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

#include <climits>
#include <cmath>    // for NAN, isfinite
#include <cstdint>  // for uint64_t, UINT64_MAX
#include <iostream>
#include <limits>   // for numeric_limits, numeric_limits<>::max_digits10
#include <set>      // for set
#include <utility>  // for pair, make_pair

CDBTTree::CDBTTree(const std::string &fname)
  : m_Filename(fname)
{
}

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
  m_FloatEntryMap[channel].insert(std::make_pair(fieldname, value));
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
  m_DoubleEntryMap[channel].insert(std::make_pair(fieldname, value));
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
  m_IntEntryMap[channel].insert(std::make_pair(fieldname, value));
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
  m_UInt64EntryMap[channel].insert(std::make_pair(fieldname, value));
}

void CDBTTree::Commit()
{
  m_Locked[MultipleEntries] = true;
  // create ntuple string
  m_TTree[MultipleEntries] = new TTree(m_TTreeName[MultipleEntries].c_str(), m_TTreeName[MultipleEntries].c_str());
  std::set<int> id_set;
  std::map<std::string, float> floatmap;
  std::map<std::string, double> doublemap;
  std::map<std::string, int> intmap;
  std::map<std::string, uint64_t> uint64map;
  intmap.insert(std::make_pair("IID", INT_MIN));
  for (auto &f_entry : m_FloatEntryMap)
  {
    id_set.insert(f_entry.first);
    for (auto &f_val : f_entry.second)
    {
      floatmap.insert(std::make_pair(f_val.first, NAN));
    }
  }
  for (auto &f_val : floatmap)
  {
    std::string fielddescriptor = f_val.first + "/F";
    m_TTree[MultipleEntries]->Branch(f_val.first.c_str(), &f_val.second, fielddescriptor.c_str());
  }

  for (auto &f_entry : m_DoubleEntryMap)
  {
    id_set.insert(f_entry.first);
    for (auto &f_val : f_entry.second)
    {
      doublemap.insert(std::make_pair(f_val.first, NAN));
    }
  }
  for (auto &f_val : doublemap)
  {
    std::string fielddescriptor = f_val.first + "/D";
    m_TTree[MultipleEntries]->Branch(f_val.first.c_str(), &f_val.second, fielddescriptor.c_str());
  }

  for (auto &i_entry : m_IntEntryMap)
  {
    id_set.insert(i_entry.first);
    for (auto &i_val : i_entry.second)
    {
      intmap.insert(std::make_pair(i_val.first, INT_MIN));
    }
  }
  for (auto &i_val : intmap)
  {
    std::string fielddescriptor = i_val.first + "/I";
    m_TTree[MultipleEntries]->Branch(i_val.first.c_str(), &i_val.second, fielddescriptor.c_str());
  }

  for (auto &i_entry : m_UInt64EntryMap)
  {
    id_set.insert(i_entry.first);
    for (auto &i_val : i_entry.second)
    {
      uint64map.insert(std::make_pair(i_val.first, UINT64_MAX));
    }
  }
  for (auto &i_val : uint64map)
  {
    std::string fielddescriptor = i_val.first + "/g";
    m_TTree[MultipleEntries]->Branch(i_val.first.c_str(), &i_val.second, fielddescriptor.c_str());
  }
  // fill ttree
  for (auto ids : id_set)
  {
    intmap["IID"] = ids;
    auto fmapiter = m_FloatEntryMap.find(ids);
    if (fmapiter != m_FloatEntryMap.end())
    {
      for (auto &f_val : fmapiter->second)
      {
        floatmap[f_val.first] = f_val.second;
      }
    }
    auto dmapiter = m_DoubleEntryMap.find(ids);
    if (dmapiter != m_DoubleEntryMap.end())
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
    for (auto &f_val : floatmap)
    {
      f_val.second = NAN;
    }
    for (auto &f_val : doublemap)
    {
      f_val.second = NAN;
    }
    for (auto &i_val : intmap)
    {
      i_val.second = INT_MIN;
    }
    for (auto &i_val : uint64map)
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
  m_TTree[SingleEntries] = new TTree(m_TTreeName[SingleEntries].c_str(), m_TTreeName[SingleEntries].c_str());
  for (auto &field : m_SingleFloatEntryMap)
  {
    std::string fielddescriptor = field.first + "/F";
    m_TTree[SingleEntries]->Branch(field.first.c_str(), &field.second, fielddescriptor.c_str());
  }
  for (auto &field : m_SingleDoubleEntryMap)
  {
    std::string fielddescriptor = field.first + "/D";
    m_TTree[SingleEntries]->Branch(field.first.c_str(), &field.second, fielddescriptor.c_str());
  }
  for (auto &field : m_SingleIntEntryMap)
  {
    std::string fielddescriptor = field.first + "/I";
    m_TTree[SingleEntries]->Branch(field.first.c_str(), &field.second, fielddescriptor.c_str());
  }
  for (auto &field : m_SingleUInt64EntryMap)
  {
    std::string fielddescriptor = field.first + "/g";
    m_TTree[SingleEntries]->Branch(field.first.c_str(), &field.second, fielddescriptor.c_str());
  }

  m_TTree[SingleEntries]->Fill();
  return;
}

void CDBTTree::Print()
{
  if (!m_FloatEntryMap.empty())
  {
    std::cout << "Number of float entries: " << m_FloatEntryMap.size() << std::endl;
    for (auto &field : m_FloatEntryMap)
    {
      std::cout << "ID: " << field.first << std::endl;
      for (auto &calibs : field.second)
      {
        std::string tmpstring = calibs.first;
        tmpstring.erase(0, 1);
        std::cout << "name " << tmpstring << " value: " << calibs.second << std::endl;
      }
    }
    std::cout << "--------------------------------------------------" << std::endl
              << std::endl;
  }
  if (!m_DoubleEntryMap.empty())
  {
    std::cout << "Number of double entries: " << m_DoubleEntryMap.size() << std::endl;
    for (auto &field : m_DoubleEntryMap)
    {
      std::cout << "ID: " << field.first << std::endl;
      for (auto &calibs : field.second)
      {
        std::string tmpstring = calibs.first;
        tmpstring.erase(0, 1);
        std::cout << "name " << tmpstring << " value: " << calibs.second << std::endl;
      }
    }
    std::cout << "--------------------------------------------------" << std::endl
              << std::endl;
  }
  if (!m_IntEntryMap.empty())
  {
    std::cout << "Number of int entries: " << m_IntEntryMap.size() << std::endl;
    for (auto &field : m_IntEntryMap)
    {
      std::cout << "ID: " << field.first << std::endl;
      for (auto &calibs : field.second)
      {
        std::string tmpstring = calibs.first;
        tmpstring.erase(0, 1);
        std::cout << "name " << tmpstring << " value: " << calibs.second << std::endl;
      }
    }
  }
  if (!m_UInt64EntryMap.empty())
  {
    std::cout << "Number of uint64 entries: " << m_UInt64EntryMap.size() << std::endl;
    for (auto &field : m_UInt64EntryMap)
    {
      std::cout << "ID: " << field.first << std::endl;
      for (auto &calibs : field.second)
      {
        std::string tmpstring = calibs.first;
        tmpstring.erase(0, 1);
        std::cout << "name " << tmpstring << " value: " << calibs.second << std::endl;
      }
    }
  }

  if (!m_SingleFloatEntryMap.empty())
  {
    std::cout << "Number of single float fields: " << m_SingleFloatEntryMap.size() << std::endl;
    for (auto &field : m_SingleFloatEntryMap)
    {
      std::string tmpstring = field.first;
      tmpstring.erase(0, 1);
      std::cout << tmpstring << " value " << field.second << std::endl;
    }
  }
  if (!m_SingleDoubleEntryMap.empty())
  {
    std::cout << "Number of single double fields: " << m_SingleDoubleEntryMap.size() << std::endl;
    // some acrobatics to restore the old state of cout after changing the precision for double printout
    std::ios oldState(nullptr);
    oldState.copyfmt(std::cout);
    std::cout.precision(std::numeric_limits<double>::max_digits10);
    for (auto &field : m_SingleDoubleEntryMap)
    {
      std::string tmpstring = field.first;
      tmpstring.erase(0, 1);
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
      tmpstring.erase(0, 1);
      std::cout << tmpstring << " value " << field.second << std::endl;
    }
  }
  if (!m_SingleUInt64EntryMap.empty())
  {
    std::cout << "Number of single uint64 fields: " << m_SingleUInt64EntryMap.size() << std::endl;
    for (auto &field : m_SingleUInt64EntryMap)
    {
      std::string tmpstring = field.first;
      tmpstring.erase(0, 1);
      std::cout << tmpstring << " value " << field.second << std::endl;
    }
  }
}

void CDBTTree::WriteCDBTTree()
{
  std::string currdir = gDirectory->GetPath();
  TFile *f = TFile::Open(m_Filename.c_str(), "RECREATE");
  for (auto ttree : m_TTree)
  {
    if (ttree != nullptr)
    {
      ttree->Write();
    }
    delete ttree;
    ttree = nullptr;
  }
  f->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

void CDBTTree::LoadCalibrations()
{
  std::string currdir = gDirectory->GetPath();

  TFile *f = TFile::Open(m_Filename.c_str());
  f->GetObject(m_TTreeName[SingleEntries].c_str(), m_TTree[SingleEntries]);
  f->GetObject(m_TTreeName[MultipleEntries].c_str(), m_TTree[MultipleEntries]);
  if (m_TTree[SingleEntries] != nullptr)
  {
    TObjArray *branches = m_TTree[SingleEntries]->GetListOfBranches();
    TIter iter(branches);
    while (TBranch *thisbranch = static_cast<TBranch *>(iter.Next()))
    {
      // this convoluted expression returns the data type of a split branch
      std::string DataType = thisbranch->GetLeaf(thisbranch->GetName())->GetTypeName();
      if (DataType == "Float_t")
      {
        auto itermap = m_SingleFloatEntryMap.insert(std::make_pair(thisbranch->GetName(), NAN));
        m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      else if (DataType == "Double_t")
      {
        auto itermap = m_SingleDoubleEntryMap.insert(std::make_pair(thisbranch->GetName(), NAN));
        m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      else if (DataType == "Int_t")
      {
        auto itermap = m_SingleIntEntryMap.insert(std::make_pair(thisbranch->GetName(), -99999));
        m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      else if (DataType == "ULong_t")
      {
        auto itermap = m_SingleUInt64EntryMap.insert(std::make_pair(thisbranch->GetName(), UINT64_MAX));
        m_TTree[SingleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
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
    m_TTree[SingleEntries]->GetEntry(0);
  }
  if (m_TTree[MultipleEntries] != nullptr)
  {
    TObjArray *branches = m_TTree[MultipleEntries]->GetListOfBranches();
    TIter iter(branches);
    std::map<std::string, float> floatvalmap;
    std::map<std::string, double> doublevalmap;
    std::map<std::string, int> intvalmap;
    std::map<std::string, uint64_t> uint64valmap;
    while (TBranch *thisbranch = static_cast<TBranch *>(iter.Next()))
    {
      // this convoluted expression returns the data type of a split branch
      std::string DataType = thisbranch->GetLeaf(thisbranch->GetName())->GetTypeName();
      if (DataType == "Float_t")
      {
        auto itermap = floatvalmap.insert(std::make_pair(thisbranch->GetName(), NAN));
        m_TTree[MultipleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      if (DataType == "Double_t")
      {
        auto itermap = doublevalmap.insert(std::make_pair(thisbranch->GetName(), NAN));
        m_TTree[MultipleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      if (DataType == "Int_t")
      {
        auto itermap = intvalmap.insert(std::make_pair(thisbranch->GetName(), INT_MIN));
        m_TTree[MultipleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
      if (DataType == "ULong_t")
      {
        auto itermap = uint64valmap.insert(std::make_pair(thisbranch->GetName(), UINT64_MAX));
        m_TTree[MultipleEntries]->SetBranchAddress(thisbranch->GetName(), &(itermap.first)->second);
      }
    }
    for (auto entry = 0; entry < m_TTree[MultipleEntries]->GetEntries(); ++entry)
    {
      for (auto &field : floatvalmap)
      {
        field.second = NAN;
      }
      for (auto &field : doublevalmap)
      {
        field.second = NAN;
      }
      for (auto &field : intvalmap)
      {
        field.second = INT_MIN;
      }
      for (auto &field : uint64valmap)
      {
        field.second = UINT64_MAX;
      }
      m_TTree[MultipleEntries]->GetEntry(entry);
      int ID = intvalmap.find("IID")->second;
      std::map<std::string, float> tmp_floatvalmap;
      for (auto &field : floatvalmap)
      {
        if (std::isfinite(field.second))
        {
          tmp_floatvalmap.insert(std::make_pair(field.first, field.second));
        }
      }
      if (!tmp_floatvalmap.empty())
      {
        m_FloatEntryMap.insert(std::make_pair(ID, tmp_floatvalmap));
      }

      std::map<std::string, double> tmp_doublevalmap;
      for (auto &field : doublevalmap)
      {
        if (std::isfinite(field.second))
        {
          tmp_doublevalmap.insert(std::make_pair(field.first, field.second));
        }
      }
      if (!tmp_doublevalmap.empty())
      {
        m_DoubleEntryMap.insert(std::make_pair(ID, tmp_doublevalmap));
      }

      std::map<std::string, int> tmp_intvalmap;
      for (auto &field : intvalmap)
      {
        if (field.second != INT_MIN && field.first != "IID")
        {
          tmp_intvalmap.insert(std::make_pair(field.first, field.second));
        }
      }
      if (!tmp_intvalmap.empty())
      {
        m_IntEntryMap.insert(std::make_pair(ID, tmp_intvalmap));
      }

      std::map<std::string, uint64_t> tmp_uint64valmap;
      for (auto &field : uint64valmap)
      {
        if (field.second != UINT64_MAX)
        {
          tmp_uint64valmap.insert(std::make_pair(field.first, field.second));
        }
      }
      if (!tmp_uint64valmap.empty())
      {
        m_UInt64EntryMap.insert(std::make_pair(ID, tmp_uint64valmap));
      }
    }
  }
  for (auto ttree : m_TTree)
  {
    delete ttree;
    ttree = nullptr;
  }
  f->Close();
  gROOT->cd(currdir.c_str());  // restore previous directory
}

float CDBTTree::GetSingleFloatValue(const std::string &name, int verbose)
{
  if (m_SingleFloatEntryMap.empty())
  {
    LoadCalibrations();
  }
  std::string fieldname = "F" + name;
  auto singleiter = m_SingleFloatEntryMap.find(fieldname);
  if (singleiter == m_SingleFloatEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " in single float calibrations" << std::endl;
      std::cout << "Existing values:" << std::endl;
      for (auto &eiter : m_SingleFloatEntryMap)
      {
        std::string tmpstring = eiter.first;
        tmpstring.erase(0, 1);
        std::cout << "name : " << tmpstring << ", value " << eiter.second
                  << std::endl;
      }
    }
    return NAN;
  }
  return singleiter->second;
}

float CDBTTree::GetFloatValue(int channel, const std::string &name, int verbose)
{
  if (m_FloatEntryMap.empty())
  {
    LoadCalibrations();
  }
  auto channelmapiter = m_FloatEntryMap.find(channel);
  if (channelmapiter == m_FloatEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find channel " << channel << " in float calibrations" << std::endl;
    }
    return NAN;
  }
  std::string fieldname = "F" + name;
  auto calibiter = channelmapiter->second.find(fieldname);
  if (calibiter == channelmapiter->second.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " among float calibrations of channel " << channel << std::endl;
    }
    return NAN;
  }
  return calibiter->second;
}

double CDBTTree::GetSingleDoubleValue(const std::string &name, int verbose)
{
  if (m_SingleDoubleEntryMap.empty())
  {
    LoadCalibrations();
  }
  std::string fieldname = "D" + name;
  auto singleiter = m_SingleDoubleEntryMap.find(fieldname);
  if (singleiter == m_SingleDoubleEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " in single double calibrations" << std::endl;
      std::cout << "Existing values:" << std::endl;
      for (auto &eiter : m_SingleDoubleEntryMap)
      {
        std::string tmpstring = eiter.first;
        tmpstring.erase(0, 1);
        std::cout << "name : " << tmpstring << ", value " << eiter.second
                  << std::endl;
      }
    }
    return NAN;
  }
  return singleiter->second;
}

double CDBTTree::GetDoubleValue(int channel, const std::string &name, int verbose)
{
  if (m_DoubleEntryMap.empty())
  {
    LoadCalibrations();
  }
  auto channelmapiter = m_DoubleEntryMap.find(channel);
  if (channelmapiter == m_DoubleEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find channel " << channel << " in double calibrations" << std::endl;
    }
    return NAN;
  }
  std::string fieldname = "D" + name;
  auto calibiter = channelmapiter->second.find(fieldname);
  if (calibiter == channelmapiter->second.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " among double calibrations for channel " << channel << std::endl;
    }
    return NAN;
  }
  return calibiter->second;
}

int CDBTTree::GetSingleIntValue(const std::string &name, int verbose)
{
  if (m_SingleIntEntryMap.empty())
  {
    LoadCalibrations();
  }
  std::string fieldname = "I" + name;
  auto singleiter = m_SingleIntEntryMap.find(fieldname);
  if (singleiter == m_SingleIntEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " in single int calibrations" << std::endl;
      std::cout << "Existing values:" << std::endl;
      for (auto &eiter : m_SingleIntEntryMap)
      {
        std::string tmpstring = eiter.first;
        tmpstring.erase(0, 1);
        std::cout << "name : " << tmpstring << ", value " << eiter.second
                  << std::endl;
      }
    }
    return INT_MIN;
  }
  return singleiter->second;
}

int CDBTTree::GetIntValue(int channel, const std::string &name, int verbose)
{
  if (m_IntEntryMap.empty())
  {
    LoadCalibrations();
  }
  auto channelmapiter = m_IntEntryMap.find(channel);
  if (channelmapiter == m_IntEntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find channel " << channel << " in int calibrations" << std::endl;
    }
    return INT_MIN;
  }
  std::string fieldname = "I" + name;
  auto calibiter = channelmapiter->second.find(fieldname);
  if (calibiter == channelmapiter->second.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " among int calibrations for channel " << channel << std::endl;
    }
    return INT_MIN;
  }
  return calibiter->second;
}

uint64_t CDBTTree::GetSingleUInt64Value(const std::string &name, int verbose)
{
  if (m_SingleUInt64EntryMap.empty())
  {
    LoadCalibrations();
  }
  std::string fieldname = "g" + name;
  auto singleiter = m_SingleUInt64EntryMap.find(fieldname);
  if (singleiter == m_SingleUInt64EntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " in single uint64 calibrations" << std::endl;
      std::cout << "Existing values:" << std::endl;
      for (auto &eiter : m_SingleUInt64EntryMap)
      {
        std::string tmpstring = eiter.first;
        tmpstring.erase(0, 1);
        std::cout << "name : " << tmpstring << ", value " << eiter.second
                  << std::endl;
      }
    }
    return UINT64_MAX;
  }
  return singleiter->second;
}

uint64_t CDBTTree::GetUInt64Value(int channel, const std::string &name, int verbose)
{
  if (m_UInt64EntryMap.empty())
  {
    LoadCalibrations();
  }
  auto channelmapiter = m_UInt64EntryMap.find(channel);
  if (channelmapiter == m_UInt64EntryMap.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find channel " << channel << " in unint64 calibrations" << std::endl;
    }
    return UINT64_MAX;
  }
  std::string fieldname = "g" + name;
  auto calibiter = channelmapiter->second.find(fieldname);
  if (calibiter == channelmapiter->second.end())
  {
    if (verbose > 0)
    {
      std::cout << "Could not find " << name << " among uint64 calibrations for channel " << channel << std::endl;
    }
    return UINT64_MAX;
  }
  return calibiter->second;
}
