#ifndef CDBOBJECTS_CDBTTREE_H
#define CDBOBJECTS_CDBTTREE_H

#include <map>
#include <set>
#include <string>

class TTree;

class CDBTTree
{
public:
  CDBTTree() = default;
  CDBTTree(const std::string &fname);
  ~CDBTTree();
  void SetFloatValue(int channel, const std::string &name, float value);
  void SetDoubleValue(int channel, const std::string &name, double value);
  void SetIntValue(int channel, const std::string &name, int value);
  void SetUInt64Value(int channel, const std::string &name, uint64_t value);
  void Commit();
  void SetSingleFloatValue(const std::string &name, float value);
  void SetSingleDoubleValue(const std::string &name, double value);
  void SetSingleIntValue(const std::string &name, int value);
  void SetSingleUInt64Value(const std::string &name, uint64_t value);
  void CommitSingle();
  void WriteCDBTTree();
  void Print();
  void LoadCalibrations();
  float GetSingleFloatValue(const std::string &name);
  float GetFloatValue(int channel, const std::string &name);
  double GetSingleDoubleValue(const std::string &name);
  double GetDoubleValue(int channel, const std::string &name);
  int GetSingleIntValue(const std::string &name);
  int GetIntValue(int channel, const std::string &name);
  uint64_t GetSingleUInt64Value(const std::string &name);
  uint64_t GetUInt64Value(int channel, const std::string &name);

private:
  enum {SingleEntries = 0, MultipleEntries = 1};
  const std::string m_TTreeName[2] = {"Single","Multiple"};
  TTree *m_TTree[2] = {nullptr};
  bool m_Locked[2] = {false};

  std::string m_Filename;
  std::set<std::string> m_EntryNameSet;
  std::map<int, std::map<std::string, float>> m_FloatEntryMap;
  std::map<std::string, float> m_SingleFloatEntryMap;
  std::map<int, std::map<std::string, double>> m_DoubleEntryMap;
  std::map<std::string, double> m_SingleDoubleEntryMap;
  std::map<int, std::map<std::string, int>> m_IntEntryMap;
  std::map<std::string, int> m_SingleIntEntryMap;
  std::map<int, std::map<std::string, uint64_t>> m_UInt64EntryMap;
  std::map<std::string, uint64_t> m_SingleUInt64EntryMap;
};

#endif
