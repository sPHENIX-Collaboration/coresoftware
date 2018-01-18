#ifndef PHParametersContainer__h
#define PHParametersContainer__h

#include <phool/PHObject.h>

#include <map>
#include <string>

class PHParameters;
class PdbParameterMapContainer;
class PHCompositeNode;

class PHParametersContainer : public PHObject
{
 public:
  typedef std::map<int, PHParameters *> Map;
  typedef Map::iterator Iterator;
  typedef Map::const_iterator ConstIterator;
  typedef std::pair<Iterator, Iterator> Range;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;

  explicit PHParametersContainer(const std::string &name = "NONE");
  virtual ~PHParametersContainer();

  void AddPHParameters(const int layer, PHParameters *params);
  const PHParameters *GetParameters(const int layer) const;
  PHParameters *GetParametersToModify(const int layer);
  int WriteToFile(const std::string &extension, const std::string &dir);
  int WriteToDB();

  void set_name(const std::string &name) { superdetectorname = name; }
  std::string Name() const { return superdetectorname; }
  //  std::pair<std::map<int, PHParameters *>::const_iterator,  std::map<int, PHParameters *>::const_iterator> GetAllParameters() {return std::make_pair(parametermap.begin(),parametermap.end());}
  ConstRange GetAllParameters() const { return std::make_pair(parametermap.begin(), parametermap.end()); }
  void Print() const;
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  int ExistDetid(const int detid) const;
  void clear() { parametermap.clear(); }
  void FillFrom(const PdbParameterMapContainer *saveparamcontainer);

 protected:
  void CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparm);
  std::string superdetectorname;
  std::map<int, PHParameters *> parametermap;
};

#endif  //PHParametersContainer__h
