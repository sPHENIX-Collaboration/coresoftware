#ifndef PHG4ParametersContainer__h
#define PHG4ParametersContainer__h

#include <phool/PHObject.h>

#include <map>
#include <string>

class PHG4Parameters;
class PdbParameterMapContainer;

class PHG4ParametersContainer: public PHObject
{
 public:
  PHG4ParametersContainer(const std::string &name = "NONE");
  virtual ~PHG4ParametersContainer();

  void AddPHG4Parameters(const int layer, PHG4Parameters *params);
  const PHG4Parameters *GetParameters(const int layer) const;
  PHG4Parameters *GetParametersToModify(const int layer);
  int WriteToFile(const std::string &extension, const std::string &dir);
  int WriteToDB();
  
  void set_name(const std::string &name) {superdetectorname = name;}
  std::string Name() const {return superdetectorname;}
  std::pair<std::map<int, PHG4Parameters *>::const_iterator,  std::map<int, PHG4Parameters *>::const_iterator> GetAllParameters() {return std::make_pair(parametermap.begin(),parametermap.end());}
  void Print() const;

 protected:
  void CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparm);
  std::string superdetectorname;
  std::map<int, PHG4Parameters *> parametermap;

};

#endif //PHG4ParametersContainer__h
