#ifndef PDBCLASSMAP_HH__
#define PDBCLASSMAP_HH__

//
// This class is a singleton wrapper around the STL-map.
// It ensures EXACTLY one existing instance of 'classmap'.
//
// Author: Matthias Messer 6/22/99

#include <cstring>
#include <map>
#include <string>

//
// The char* comparison operator, required by the STL map
//
template <typename T> 
struct strless : 
  std::binary_function<T, T, bool> {
  bool operator()(const T& x, const T& y) const
  {
    return strcmp(x.c_str(),y.c_str()) < 0;
  }
};

template <typename T> 
class PdbClassMap 
{
public:

  static PdbClassMap *instance();
  virtual ~PdbClassMap();

  T*& operator [] (const std::string & className) { return _map[className]; }
  typename std::map<const std::string, T*, strless<const std::string> >::iterator find(const std::string & className) { return _map.find(className); }
  typename std::map<const std::string, T*, strless<const std::string> >::iterator end() { return _map.end(); }
  void erase(const std::string & className);

protected:
  PdbClassMap();
  
private:
  static PdbClassMap * _instance;
  std::map<const std::string, T*, strless<const std::string> > _map;
};

template <typename T> 
PdbClassMap<T>::PdbClassMap()
{
}

template <typename T> 
void PdbClassMap<T>::erase(const std::string & className)
{
  _map.erase(_map.find(className));
  // if we have no more entries - delete ourselves
  if (_map.size() == 0)
    {
      delete _instance;
      _instance = NULL;
    }
}

template <typename T> 
PdbClassMap<T>::~PdbClassMap()
{
  while( _map.begin() != _map.end())
    {
      delete _map.begin()->second;
      _map.erase(_map.begin());
    }
}

template <typename T> 
PdbClassMap<T>* PdbClassMap<T>::instance()
{
  if (!_instance)
    {
      _instance = new PdbClassMap<T>();
    }
  return _instance;
}

template <typename T> PdbClassMap<T>* PdbClassMap<T>::_instance = NULL;

#endif /* __PDBCLASSMAP_HH__ */
