// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHPARAMETER_PHPARAMETERSCONTAINER_H
#define PHPARAMETER_PHPARAMETERSCONTAINER_H

#include <phool/PHObject.h>

#include <map>
#include <string>
#include <utility>

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
  ~PHParametersContainer() override;

  void AddPHParameters(const int detid, PHParameters *params);
  const PHParameters *GetParameters(const int detid) const;
  PHParameters *GetParametersToModify(const int detid);
  int WriteToFile(const std::string &extension, const std::string &dir);
  int WriteToDB();

  void set_name(const std::string &name) { superdetectorname = name; }
  std::string Name() const { return superdetectorname; }
  //  std::pair<std::map<int, PHParameters *>::const_iterator,  std::map<int, PHParameters *>::const_iterator> GetAllParameters() {return std::make_pair(parametermap.begin(),parametermap.end());}
  ConstRange GetAllParameters() const { return std::make_pair(parametermap.begin(), parametermap.end()); }
  void Print(Option_t *option = "") const override;
  void SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  void UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename);
  int ExistDetid(const int detid) const;
  void clear() { parametermap.clear(); }
  void FillFrom(const PdbParameterMapContainer *saveparamcontainer);
  void CreateAndFillFrom(const PdbParameterMapContainer *saveparamcontainer, const std::string &name);

 private:
  void CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparm);
  void UpdatePdbParameterMapContainer(PdbParameterMapContainer *myparm);
  std::string superdetectorname;
  std::map<int, PHParameters *> parametermap;
};

#endif  //  PHPARAMETER_PHPARAMETERSCONTAINER_H
