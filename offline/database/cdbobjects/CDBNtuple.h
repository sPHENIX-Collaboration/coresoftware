#ifndef CDBOBJECTS_CDBNTUPLE_H
#define CDBOBJECTS_CDBNTUPLE_H

#include <map>
#include <set>
#include <string>

class TNtuple;

class CDBNtuple
{
public:
  CDBNtuple() = default;
  CDBNtuple(const std::string &fname);
  ~CDBNtuple();
  void SetValue(int channel, const std::string &name, float value);
  void Commit();
  void SetSingleValue(const std::string &name, float value);
  void CommitSingle();
  void WriteCDBNtuple();
  void Print();
  void LoadCalibrations();
  float GetSingleValue(const std::string &name);
  float GetValue(int channel, const std::string &name);

private:
  enum {SingleEntries = 0, MultipleEntries = 1};
  const std::string m_TNtupleName[2] = {"Single","Multiple"};
  TNtuple *m_TNtuple[2] = {nullptr};
  bool m_Locked[2] = {false};

  std::string m_Filename;
  std::set<std::string> m_EntryNameSet;
  std::map<int, std::map<std::string, float>> m_EntryMap;
  std::map<std::string, float> m_SingleEntryMap;
};

#endif
