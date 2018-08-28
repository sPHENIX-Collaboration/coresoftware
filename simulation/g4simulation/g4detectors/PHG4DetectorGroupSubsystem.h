#ifndef PHG4DetectorGroupSubsystem_h
#define PHG4DetectorGroupSubsystem_h

#include <g4main/PHG4Subsystem.h>

#include <map>
#include <set>
#include <string>

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

  virtual ~PHG4DetectorGroupSubsystem() {}
// stupid rootcint does not support final keyword
#ifndef __CINT__
  int Init(PHCompositeNode *) final;
  int InitRun(PHCompositeNode *) final;
#else
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
#endif

  virtual int InitRunSubsystem(PHCompositeNode *)
  {
    return 0;
  }
  virtual int InitSubsystem(PHCompositeNode *) { return 0; }
  void OverlapCheck(const bool chk = true) { overlapcheck = chk; }
  bool CheckOverlap() const { return overlapcheck; }
  PHParametersContainer *GetParamsContainer() const { return paramscontainer; }
  // Get/Set parameters from macro
  void set_double_param(const int detid, const std::string &name, const double dval);
  double get_double_param(const int detid, const std::string &name) const;
  void set_int_param(const int detid, const std::string &name, const int ival);
  int get_int_param(const int detid, const std::string &name) const;
  void set_string_param(const int detid, const std::string &name, const std::string &sval);
  std::string get_string_param(const int detid, const std::string &name) const;

  void UseDB(const int i = 1) { usedb = i; }
  int ReadDB() const { return usedb; }
  FILE_TYPE get_filetype() const { return filetype; }
  void UseCalibFiles(const FILE_TYPE ftyp) { filetype = ftyp; }
  int SaveParamsToDB();
  int ReadParamsFromDB(const std::string &name, const int issuper);
  int SaveParamsToFile(const FILE_TYPE ftyp);
  int ReadParamsFromFile(const std::string &name, const FILE_TYPE ftyp, const int issuper);
  void SetCalibrationFileDir(const std::string &calibdir) { calibfiledir = calibdir; }
  void UpdateParametersWithMacro();

  void SetActive(const int detid, const int i);
  void SetActive(const int i = 1);
  void SetAbsorberActive(const int detid, const int i);
  void SetAbsorberActive(const int i = 1);
  void SetAbsorberTruth(const int detid, const int i);
  void SetAbsorberTruth(const int i = 1);
  void BlackHole(const int detid, const int i);
  void BlackHole(const int i = 1);
  void SuperDetector(const std::string &name);
  const std::string SuperDetector() const { return superdetector; }
  int GetLayer() const { return layer; }
  virtual void SetDefaultParameters() = 0;  // this one has to be implemented by the daughter
 protected:                                 // those cannot be executed on the cmd line
  PHG4DetectorGroupSubsystem(const std::string &name = "GenericSubsystem", const int lyr = 0);
  // these initialize the defaults and add new entries to the
  // list of variables. This should not be possible from the macro to
  // prevent abuse (this makes the list of possible parameters deterministic)
  void InitializeParameters();
  void AddDetId(const int i) { layers.insert(i); }
  std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> GetDetIds() const
  {
    return std::make_pair(layers.begin(), layers.end());
  }
  void set_default_double_param(const int detid, const std::string &name, const double dval);
  void set_default_int_param(const int detid, const std::string &name, const int ival);
  void set_default_string_param(const int detid, const std::string &name, const std::string &sval);
  int BeginRunExecuted() const { return beginrunexecuted; }
  void PrintDefaultParams() const;
  void PrintMacroParams() const;

 private:
  PHParametersContainer *paramscontainer;
  PHParametersContainer *paramscontainer_default;
  PHCompositeNode *savetopNode;
  bool overlapcheck;
  int layer;
  int usedb;
  int beginrunexecuted;
  FILE_TYPE filetype;
  std::string superdetector;
  std::string calibfiledir;

  std::set<int> layers;

  std::map<int, std::map<const std::string, double>> dparams;
  std::map<int, std::map<const std::string, int>> iparams;
  std::map<int, std::map<const std::string, std::string>> cparams;

  std::map<int, std::map<const std::string, double>> default_double;
  std::map<int, std::map<const std::string, int>> default_int;
  std::map<int, std::map<const std::string, std::string>> default_string;
};

#endif
