#ifndef PHG4ParametersContainer__h
#define PHG4ParametersContainer__h

#include <phool/PHObject.h>

#include <map>
#include <string>

class PHG4Parameters;
class PdbParameterMapContainer;
class PHCompositeNode;

//! This class is deprecated, please use PHG4Parameter instead.
//! See https://github.com/sPHENIX-Collaboration/coresoftware/pull/405
class PHG4ParametersContainer : public PHObject
{
 public:
  typedef std::map<int, PHG4Parameters *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  explicit PHG4ParametersContainer(const std::string &name = "NONE");
  virtual ~PHG4ParametersContainer();

  void AddPHG4Parameters(const int layer, PHG4Parameters *params);
  const PHG4Parameters *GetParameters(const int layer) const;
  PHG4Parameters *GetParametersToModify(const int layer);
  int WriteToFile(const std::string &extension, const std::string &dir);
  int WriteToDB();

  void set_name(const std::string &name) { superdetectorname = name; }
  std::string Name() const { return superdetectorname; }
  //  std::pair<std::map<int, PHG4Parameters *>::const_iterator,  std::map<int, PHG4Parameters *>::const_iterator> GetAllParameters() {return std::make_pair(parametermap.begin(),parametermap.end());}
  ConstRange GetAllParameters() const { return std::make_pair(parametermap.begin(), parametermap.end()); }
  void Print() const;
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  int ExistDetid(const int detid) const;
  void clear() { parametermap.clear(); }
  void FillFrom(const PdbParameterMapContainer *saveparamcontainer);
  void CreateAndFillFrom(const PdbParameterMapContainer *saveparamcontainer, const std::string &name);

 protected:
  void CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparm);
  std::string superdetectorname;
  std::map<int, PHG4Parameters *> parametermap;
};

#endif  //PHG4ParametersContainer__h
