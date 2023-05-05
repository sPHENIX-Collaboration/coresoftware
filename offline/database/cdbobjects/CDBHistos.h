#ifndef CDBOBJECTS_CDBHISTOS_H
#define CDBOBJECTS_CDBHISTOS_H

#include <map>
#include <string>

class TH1;

class CDBHistos
{
 public:
  CDBHistos() = default;
  explicit CDBHistos(const std::string &fname);
  ~CDBHistos();
  void WriteCDBHistos();
  void Print() const;
  void LoadCalibrations();
  void registerHisto(TH1 *h1);
  TH1 *getHisto(const std::string &name);

 private:
  std::string m_Filename;
  std::map<std::string, TH1 *> m_HistoMap;
};

#endif
