
#ifndef __EVALLINKSV1_H__
#define __EVALLINKSV1_H__

#include "EvalLinks.h"

#include <iostream>
#include <string>
#include <set>
#include <map>

class EvalLinksV1 : public EvalLinks {

public:

  EvalLinksV1(const std::string& left_name = "left",
	      const std::string& right_name = "right",
	      const std::string& weight_name = "weight");
  virtual ~EvalLinksV1() {clear();}

  void identify(std::ostream& os = std::cout) const;
  void Reset() {clear();}
  int isValid() const {return 1;}
  
  // modifiers
  void set_names(const std::string &left_name,
		 const std::string &right_name,
		 const std::string &weight_name);
  void link(unsigned long long left_id, unsigned long long right_id, double weight);
  void unlink(unsigned long long left_id, unsigned long long right_id);
  void unlink_subleading();
  void clear();

  // status
  size_t size() const;
  std::string get_name_left() const {return _left_name;}
  std::string get_name_right() const {return _right_name;}
  std::string get_name_weight() const {return _weight_name;} 
  bool        has_link(unsigned long long left_id, unsigned long long right_id) const;
  double  get_weight(unsigned long long left_id, unsigned long long right_id) const;

  // status for missing links
  bool               is_null_id(unsigned long long id) const {return (id == NULLID);}
  unsigned long long get_null_id() const {return NULLID;}
  
  // access all associations
  std::set<unsigned long long> left(unsigned long long right_id) const;
  std::set<unsigned long long> right(unsigned long long left_id) const;

  // access best weight association
  unsigned long long max_left(unsigned long long right_id) const;
  unsigned long long max_right(unsigned long long left_id) const;

private:

  static const unsigned long long NULLID;
  
  bool stale() const {return _stale;}
  void refresh() const;

  void calc_max_left(unsigned long long right_id) const;
  void calc_max_right(unsigned long long left_id) const;

  std::string _left_name;  //< left object container names (e.g. SvtxTrackMap)
  std::string _right_name; //< right object containter names (e.g. G4TrughtInfo)
  std::string _weight_name; //< connection weight meaning (e.g. nhits)
  
  /// storage for (left id,right id) => weight value
  std::map<std::pair<unsigned long long,unsigned long long>, double> _links;

#ifndef __CINT____ // hide from dictionary generation
  mutable bool _stale; //! exclude from ROOT I/O
  mutable std::multimap<unsigned long long,unsigned long long> _left_right_mmap; //! exclude from ROOT I/O
  mutable std::multimap<unsigned long long,unsigned long long> _right_left_mmap; //! exclude from ROOT I/O
  mutable std::map<unsigned long long,unsigned long long> _left_right_map; //! exclude from ROOT I/O
  mutable std::map<unsigned long long,unsigned long long> _right_left_map; //! exclude from ROOT I/O
#endif // __CINT__

  ClassDef(EvalLinksV1,1)
};

#endif // __EVALLINKSV1_H__
