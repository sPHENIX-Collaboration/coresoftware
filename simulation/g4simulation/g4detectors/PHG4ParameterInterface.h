#ifndef PHG4ParameterInterface__H
#define PHG4ParameterInterface__H

#include <map>
#include <string>

class PHCompositeNode;
class PHG4Parameters;

//! This class is deprecated, please use PHG4Parameter instead.
//! See https://github.com/sPHENIX-Collaboration/coresoftware/pull/405
class PHG4ParameterInterface
{
 public:
  PHG4ParameterInterface(const std::string &name);
  virtual ~PHG4ParameterInterface(){}

  void set_paramname(const std::string &name);
  virtual void  SetDefaultParameters() = 0;

 // Get/Set parameters from macro
  void set_double_param(const std::string &name, const double dval);
  double get_double_param(const std::string &name) const;
  void set_int_param(const std::string &name, const int ival);
  int get_int_param(const std::string &name) const;
  void set_string_param(const std::string &name, const std::string &sval);
  std::string get_string_param(const std::string &name) const;

  void UpdateParametersWithMacro();
  void SaveToNodeTree(PHCompositeNode *runNode, const std::string &nodename);
  void PutOnParNode(PHCompositeNode *parNode, const std::string &nodename);
 protected:
  void set_default_double_param( const std::string &name, const double dval);
  void set_default_int_param( const std::string &name, const int ival);
  void set_default_string_param( const std::string &name, const std::string &sval);
  void InitializeParameters();

 private:
  PHG4Parameters *params;
  std::map<const std::string, double> dparams;
  std::map<const std::string, int> iparams;
  std::map<const std::string, std::string> cparams;

  std::map<const std::string, double> default_double;
  std::map<const std::string, int> default_int;
  std::map<const std::string, std::string> default_string;

};

#endif
