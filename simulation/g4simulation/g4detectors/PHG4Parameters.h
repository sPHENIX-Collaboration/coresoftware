#ifndef PHG4PARAMETERS_H
#define PHG4PARAMETERS_H

#include <phool/PHObject.h>

#include <map>
#include <string>

class PdbParameterMap;
class PHCompositeNode;

// contains parameters in our units,
// convert to G4 units inside get access methods

class PHG4Parameters: public PHObject
{
 public:

  typedef std::map<const std::string, double> dMap;
  typedef std::map<const std::string, int> iMap;
  typedef std::map<const std::string, std::string> strMap;

 PHG4Parameters(const std::string &name): pdbparam(NULL),detname(name) {}
  virtual ~PHG4Parameters() {}

  void print() const;

  void set_int_param(const std::string &name, const int ival);
  int get_int_param(const std::string &name) const;

  void set_double_param(const std::string &name, const double dval);
  double get_double_param(const std::string &name) const;

  void set_string_param(const std::string &name, const std::string &str);
  std::string get_string_param(const std::string &name) const;

  void set_name(const std::string &name) {detname = name;}

  void FillFrom(const PdbParameterMap *saveparams);
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  void WriteToDB();

 protected:
  void printint() const;
  void printdouble() const;
  void printstring() const;
  PdbParameterMap *pdbparam;
  std::string detname;
  std::map<const std::string, double> doubleparams;
  std::map<const std::string, int> intparams;
  std::map<const std::string, std::string> stringparams;

};

#endif
