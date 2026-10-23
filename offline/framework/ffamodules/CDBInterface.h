// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_CDBINTERFACE_H
#define FFAMODULES_CDBINTERFACE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>  // for uint64_t
#include <map>
#include <set>
#include <string>
#include <tuple>  // for tuple

class SphenixClient;

class CDBInterface : public SubsysReco
{
 public:
  static CDBInterface *instance();

  ~CDBInterface() override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  int UpdateRunNode(PHCompositeNode *topNode) override;

  void Disable() { disable = true; }
  void Enable() { disable = false; }

  void Disable_default() { disable_default = true; }
  void Enable_default() { disable_default = false; }

  std::string getUrl(const std::string &domain, const std::string &filename = "");

  void DumpCalibrations(const std::string &filename);
  void ReadCalibrationsFromFile(const std::string &filename);

  //added for in-run transmission of cdb information
  // returns a copy (not a lot of additional memory and ran right before CDBInterface::End()
  // rather this than risk having a reference around
  
  std::set<std::tuple<std::string, std::string, uint64_t>> getUrls() const { return m_UrlVector; }

  // safety check for other code interacting with CDBInterface
  static bool exists() { return __instance != nullptr; }

private:
  CDBInterface(const std::string &name = "CDBInterface");

  static CDBInterface *__instance;
  SphenixClient *cdbclient{nullptr};
  bool disable{false};
  bool disable_default{false};
  bool m_Read_From_File_Flag{false};
  std::map<std::string, std::string> m_Payload_Url_Cache;
  std::set<std::tuple<std::string, std::string, uint64_t>> m_UrlVector;
};

#endif  // FFAMODULES_CDBINTERFACE_H
