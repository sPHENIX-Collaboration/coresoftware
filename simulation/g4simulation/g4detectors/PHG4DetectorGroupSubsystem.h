// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4DETECTORGROUPSUBSYSTEM_H
#define G4DETECTORS_PHG4DETECTORGROUPSUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <map>
#include <set>
#include <string>
#include <utility>  // for make_pair, pair

class PHCompositeNode;
class PHParametersContainer;

class PHG4DetectorGroupSubsystem : public PHG4Subsystem
{
 public:
  enum FILE_TYPE
  {
    none = 0,
    xml = 1,
    root = 2
  };

  ~PHG4DetectorGroupSubsystem() override {}
  int Init(PHCompositeNode *) final;
  int InitRun(PHCompositeNode *) final;

  virtual int InitRunSubsystem(PHCompositeNode *)
  {
    return 0;
  }
  virtual int InitSubsystem(PHCompositeNode *) { return 0; }
  void OverlapCheck(const bool chk = true) { m_OverlapCheckFlag = chk; }
  bool CheckOverlap() const { return m_OverlapCheckFlag; }
  PHParametersContainer *GetParamsContainer() const { return m_ParamsContainer; }
  // Get/Set parameters from macro
  void set_double_param(const int detid, const std::string &name, const double dval);
  double get_double_param(const int detid, const std::string &name) const;
  void set_int_param(const int detid, const std::string &name, const int ival);
  int get_int_param(const int detid, const std::string &name) const;
  void set_string_param(const int detid, const std::string &name, const std::string &sval);
  std::string get_string_param(const int detid, const std::string &name) const;

  void UseDB(const int i = 1) { m_UseDBFlag = i; }
  int ReadDB() const { return m_UseDBFlag; }
  FILE_TYPE get_filetype() const { return m_FileType; }
  void UseCalibFiles(const FILE_TYPE ftyp) { m_FileType = ftyp; }
  int SaveParamsToDB();
  int ReadParamsFromDB(const std::string &name, const int issuper);
  int SaveParamsToFile(const FILE_TYPE ftyp);
  int ReadParamsFromFile(const std::string &name, const FILE_TYPE ftyp, const int issuper);
  void SetCalibrationFileDir(const std::string &calibdir) { m_CalibFileDir = calibdir; }
  void UpdateParametersWithMacro();

  void SetActive(const int detid, const int i);
  void SetActive(const int i = 1);
  void SetAbsorberActive(const int detid, const int i);
  void SetAbsorberActive(const int i = 1);
  void SetAbsorberTruth(const int detid, const int i);
  void SetAbsorberTruth(const int i = 1);
  void SetSupportActive(const int detid, const int i = 1);
  void SetSupportActive(const int i = 1);

  void BlackHole(const int detid, const int i);
  void BlackHole(const int i = 1);
  void SuperDetector(const std::string &name);
  const std::string SuperDetector() const { return m_SuperDetector; }
  int GetLayer() const { return m_Layer; }
  virtual void SetDefaultParameters() = 0;  // this one has to be implemented by the daughter

 protected:  // those cannot be executed on the cmd line
  PHG4DetectorGroupSubsystem(const std::string &name = "GenericSubsystem", const int lyr = 0);
  // these initialize the defaults and add new entries to the
  // list of variables. This should not be possible from the macro to
  // prevent abuse (this makes the list of possible parameters deterministic)
  void InitializeParameters();
  void AddDetId(const int i) { m_LayerSet.insert(i); }
  std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> GetDetIds() const
  {
    return std::make_pair(m_LayerSet.begin(), m_LayerSet.end());
  }
  void set_default_double_param(const int detid, const std::string &name, const double dval);
  void set_default_int_param(const int detid, const std::string &name, const int ival);
  void set_default_string_param(const int detid, const std::string &name, const std::string &sval);
  int BeginRunExecuted() const { return m_BeginRunExecutedFlag; }
  void PrintDefaultParams() const;
  void PrintMacroParams() const;

 private:
  PHParametersContainer *m_ParamsContainer = nullptr;
  PHParametersContainer *m_ParamsContainerDefault = nullptr;
  PHCompositeNode *m_SaveTopNode = nullptr;
  bool m_OverlapCheckFlag = false;
  int m_Layer = 0;
  int m_UseDBFlag = 0;
  int m_BeginRunExecutedFlag = 0;
  FILE_TYPE m_FileType = PHG4DetectorGroupSubsystem::none;
  std::string m_SuperDetector = "NONE";
  std::string m_CalibFileDir = "./";

  std::set<int> m_LayerSet;

  std::map<int, std::map<const std::string, double>> m_MacroDoubleParamsMap;
  std::map<int, std::map<const std::string, int>> m_MacroIntegerParamsMap;
  std::map<int, std::map<const std::string, std::string>> m_MacroStringParamsMap;

  std::map<int, std::map<const std::string, double>> m_DefaultDoubleParamsMap;
  std::map<int, std::map<const std::string, int>> m_DefaultIntegerParamsMap;
  std::map<int, std::map<const std::string, std::string>> m_DefaultStringParamsMap;
};

#endif
