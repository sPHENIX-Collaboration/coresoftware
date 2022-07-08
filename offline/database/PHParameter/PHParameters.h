// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHPARAMETER_PHPARAMETERS_H
#define PHPARAMETER_PHPARAMETERS_H

#include <phool/PHObject.h>

#include <cstddef>
#include <map>
#include <string>
#include <utility>

class PdbParameterMap;
class PdbParameterMapContainer;
class PHCompositeNode;

// contains parameters in our units,
// convert to G4 units inside get access methods

class PHParameters : public PHObject
{
 public:
  typedef std::map<const std::string, double> dMap;
  typedef dMap::const_iterator dIter;
  typedef std::map<const std::string, int> iMap;
  typedef iMap::const_iterator iIter;
  typedef std::map<const std::string, std::string> strMap;
  typedef std::map<const std::string, std::string>::const_iterator strIter;

  explicit PHParameters(const std::string &name)
    : m_Detector(name)
  {
  }

  PHParameters(const PHParameters &params, const std::string &name);

  ~PHParameters() override;

  void Reset() override;

  void Print(Option_t *option = "") const override;

  //! hash of binary information for checking purpose
  size_t get_hash() const;

  void set_int_param(const std::string &name, const int ival);
  int get_int_param(const std::string &name) const;
  bool exist_int_param(const std::string &name) const;
  std::pair<std::map<const std::string, int>::const_iterator, std::map<const std::string, int>::const_iterator> get_all_int_params() { return std::make_pair(m_IntParMap.begin(), m_IntParMap.end()); }

  void set_double_param(const std::string &name, const double dval);
  double get_double_param(const std::string &name) const;
  bool exist_double_param(const std::string &name) const;
  std::pair<std::map<const std::string, double>::const_iterator, std::map<const std::string, double>::const_iterator> get_all_double_params() { return std::make_pair(m_DoubleParMap.begin(), m_DoubleParMap.end()); }

  void set_string_param(const std::string &name, const std::string &str);
  std::string get_string_param(const std::string &name) const;
  bool exist_string_param(const std::string &name) const;
  std::pair<std::map<const std::string, std::string>::const_iterator, std::map<const std::string, std::string>::const_iterator> get_all_string_params() { return std::make_pair(m_StringParMap.begin(), m_StringParMap.end()); }

  void set_name(const std::string &name) { m_Detector = name; }
  std::string Name() const { return m_Detector; }

  void FillFrom(const PdbParameterMap *saveparams);
  void FillFrom(const PdbParameterMapContainer *saveparamcontainer, const int layer);
  void FillFrom(const PHParameters *saveparams);
  // save parameters on node tree
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  // save parameters in container on node tree
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename, const int layer);

  // update parameters on node tree (in case the subsystem modified them)
  void UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  void UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename, const int layer);

  int WriteToDB();
  int ReadFromDB();
  int ReadFromCDB(const std::string &domain);
  int ReadFromDB(const std::string &name, const int layer);
  int WriteToCDBFile(const std::string &filename);
  int WriteToFile(const std::string &extension, const std::string &dir = ".");

  //! simple read without super detector and layer structures
  inline int ReadFromFile(const std::string &name, const std::string &extension, const std::string &dir = ".")
  {
    return ReadFromFile(name, extension, 0, 0, dir);
  }
  int ReadFromCDBFile(const std::string &name);
  //! Fully fledged read
  int ReadFromFile(const std::string &name, const std::string &extension, const int layer, const int issuper, const std::string &dir = ".");
  void CopyToPdbParameterMap(PdbParameterMap *myparm);

  void printint() const;
  void printdouble() const;
  void printstring() const;

 private:
  unsigned int ConvertStringToUint(const std::string &str) const;
  std::string m_Detector;
  dMap m_DoubleParMap;
  iMap m_IntParMap;
  strMap m_StringParMap;

  //No Class Def since this class is not intended to be persistent
};

#endif  // PHPARAMETER_PHPARAMETERS_H
