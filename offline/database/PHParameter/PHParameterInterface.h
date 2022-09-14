// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHPARAMETER_PHPARAMETERINTERFACE_H
#define PHPARAMETER_PHPARAMETERINTERFACE_H

#include <map>
#include <string>

class PHCompositeNode;
class PHParameters;

class PHParameterInterface
{
 public:
  PHParameterInterface(const std::string &name);
  // PHParameterInterface contains pointer to memory
  // copy ctor and = operator need explicit implementation, do just delete it here
  PHParameterInterface(const PHParameterInterface &) = delete;
  PHParameterInterface &operator=(PHParameterInterface const &) = delete;

  virtual ~PHParameterInterface();

  void set_paramname(const std::string &name);
  virtual void SetDefaultParameters() = 0;

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
  void Print() const;

 protected:
  void set_default_double_param(const std::string &name, const double dval);
  void set_default_int_param(const std::string &name, const int ival);
  void set_default_string_param(const std::string &name, const std::string &sval);
  void InitializeParameters();

 private:
  PHParameters *m_Params = nullptr;

  bool m_Locked = false;
  std::map<const std::string, double> m_DoubleParMap;
  std::map<const std::string, int> m_IntParMap;
  std::map<const std::string, std::string> m_StringParMap;

  std::map<const std::string, double> m_DefaultDoubleParMap;
  std::map<const std::string, int> m_DefaultIntParMap;
  std::map<const std::string, std::string> m_DefaultStringParMap;
};

#endif  // PHPARAMETER_PHPARAMETERINTERFACE_H
