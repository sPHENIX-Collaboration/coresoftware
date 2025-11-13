#ifndef CDBOBJECTS_CDBTTREE_H
#define CDBOBJECTS_CDBTTREE_H

#include <cstdint>
#include <map>
#include <string>

class TTree;

class CDBTTree
{
 public:
  CDBTTree() = default;
  explicit CDBTTree(const std::string &fname);
  ~CDBTTree();
  static void SetVerbosity(int v) { verbosity = v; }
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
  void WriteSingleCDBTTree();
  void WriteMultipleCDBTTree();
  void Print();
  void SetFilename(const std::string &fname) { m_Filename = fname; }
  void LoadCalibrations();
  float GetSingleFloatValue(const std::string &name, int verbose = 0);
  float GetFloatValue(int channel, const std::string &name, int verbose = 0);
  double GetSingleDoubleValue(const std::string &name, int verbose = 0);
  double GetDoubleValue(int channel, const std::string &name, int verbose = 0);
  int GetSingleIntValue(const std::string &name, int verbose = 0);
  int GetIntValue(int channel, const std::string &name, int verbose = 0);
  uint64_t GetSingleUInt64Value(const std::string &name, int verbose = 0);
  uint64_t GetUInt64Value(int channel, const std::string &name, int verbose = 0);

  const auto &GetFloatEntryMap() const { return m_FloatEntryMap; }
  const auto &GetDoubleEntryMap() const { return m_DoubleEntryMap; }
  const auto &GetIntEntryMap() const { return m_IntEntryMap; }
  const auto &GetUInt64EntryMap() const { return m_UInt64EntryMap; }

  const auto &GetSingleFloatEntryMap() const { return m_SingleFloatEntryMap; }
  const auto &GetSingleDoubleEntryMap() const { return m_SingleDoubleEntryMap; }
  const auto &GetSingleIntEntryMap() const { return m_SingleIntEntryMap; }
  const auto &GetSingleUInt64EntryMap() const { return m_SingleUInt64EntryMap; }

 private:
  enum
  {
    SingleEntries = 0,
    MultipleEntries = 1
  };
  const std::string m_TTreeName[2] = {"Single", "Multiple"};
  TTree *m_TTree[2] = {nullptr};
  static int verbosity;
  bool m_Locked[2] = {false};

  std::string m_Filename;
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
