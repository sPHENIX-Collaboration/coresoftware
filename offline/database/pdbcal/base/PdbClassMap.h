#ifndef PDBCAL_BASE_PDBCLASSMAP_H
#define PDBCAL_BASE_PDBCLASSMAP_H

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
    return strcmp(x,y) < 0;
  }
};

template <typename T> 
class PdbClassMap 
{
public:

  static PdbClassMap *instance();
  virtual ~PdbClassMap();

  T*& operator [] (const char * className) { return _map[className]; }
  typename std::map<const char*, T*, strless<const char*> >::iterator find(const char * className) { return _map.find(className); }
  typename std::map<const char*, T*, strless<const char*> >::iterator end() { return _map.end(); }
  void erase(const char * className);

protected:
  PdbClassMap();
  
private:
  static PdbClassMap * _instance;
  std::map<const char*, T*, strless<const char*> > _map;
};

template <typename T> 
PdbClassMap<T>::PdbClassMap()
{
}

template <typename T> 
void PdbClassMap<T>::erase(const char * className)
{
  _map.erase(_map.find(className));
  // if we have no more entries - delete ourselves
  if (_map.size() == 0)
    {
      delete _instance;
      _instance = 0;
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

template <typename T> PdbClassMap<T>* PdbClassMap<T>::_instance = 0;

#endif /* PDBCAL_BASE_PDBCLASSMAP_H */
