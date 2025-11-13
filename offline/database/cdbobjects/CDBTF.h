#ifndef CDBOBJECTS_CDBTF_H
#define CDBOBJECTS_CDBTF_H

#include <map>
#include <string>

class TF1;

class CDBTF
{
 public:
  CDBTF() = default;
  explicit CDBTF(const std::string &fname);
  ~CDBTF();
  void WriteCDBTF();
  void Print() const;
  void LoadCalibrations();
  void registerTF(TF1 *t1);
  TF1 *getTF(const std::string &name, bool printerror = true);
  const auto &GetTFMap() const {return m_TFMap;}

 private:
  std::string m_Filename;
  std::map<std::string, TF1 *> m_TFMap;
};

#endif
