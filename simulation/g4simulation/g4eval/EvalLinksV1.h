
#ifndef __EVALLINKSV1_H__
#define __EVALLINKSV1_H__

#include "EvalLinks.h"

#include <string>
#include <set>
#include <map>

class EvalLinksV1 : public EvalLinks {

public:

  EvalLinksV1(std::string left_name, std::string right_name);
  virtual ~EvalLinksV1() {clear();}

  void Reset() {clear();}
  
  // modifiers
  void set_names(std::string left_name, std::string right_name);
  void link(unsigned int left_id, unsigned int right_id, float purity);
  void unlink(unsigned int left_id, unsigned int right_id);
  void clear();

  // status
  void   print() const;
  size_t size() const;
  bool   has_link(unsigned int left_id, unsigned int right_id) const;
  float  get_purity(unsigned int left_id, unsigned int right_id) const;
  
  // access all associations
  std::set<unsigned int> left(unsigned int right_id) const;
  std::set<unsigned int> right(unsigned int left_id) const;

  // access best purity association
  unsigned int max_left(unsigned int right_id) const {return _right_left_map[right_id];}
  unsigned int max_right(unsigned int left_id) const {return _left_right_map[left_id];}

private:

  bool stale() const {return _stale;}
  void refresh() const;

  void calc_max_left(unsigned int right_id) const;
  void calc_max_right(unsigned int left_id) const;

  std::string _left_name;  //< left object container names (e.g. SvtxTrackMap)
  std::string _right_name; //< right object containter names (e.g. PHG4ParticleMap)
  
  /// storage for (left id,right id) => purity value
  std::map<std::pair<unsigned int,unsigned int>, float> _links;

#ifndef __CINT____ // hide from dictionary generation
  mutable bool _stale; //! exclude from ROOT I/O
  mutable std::multimap<unsigned int,unsigned int> _left_right_mmap; //! exclude from ROOT I/O
  mutable std::multimap<unsigned int,unsigned int> _right_left_mmap; //! exclude from ROOT I/O
  mutable std::map<unsigned int,unsigned int> _left_right_map; //! exclude from ROOT I/O
  mutable std::map<unsigned int,unsigned int> _right_left_map; //! exclude from ROOT I/O
#endif // __CINT__

  ClassDef(EvalLinksV1,1)
};

#endif // __EVALLINKSV1_H__
