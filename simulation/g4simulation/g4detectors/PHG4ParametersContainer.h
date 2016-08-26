#ifndef PHG4ParametersContainer__h
#define PHG4ParametersContainer__h

#include <phool/PHObject.h>

#include <map>
#include <string>

class PHG4Parameters;

class PHG4ParametersContainer: public PHObject
{
 public:
  PHG4ParametersContainer() {}
  virtual ~PHG4ParametersContainer();

  void AddPHG4Parameters(const int layer, PHG4Parameters *params);
  const PHG4Parameters *GetParameters(const int layer) const;
  PHG4Parameters *GetParametersToModify(const int layer);

  void set_name(const std::string &name) {superdetectorname = name;}
  std::string Name() const {return superdetectorname;}

 protected:
  std::string superdetectorname;
  std::map<int, PHG4Parameters *> parametermap;

};

#endif //PHG4ParametersContainer__h
