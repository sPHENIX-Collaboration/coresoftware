
#ifndef __EVALLINKS_H__
#define __EVALLINKS_H__

#include <phool/PHObject.h>

#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include <limits.h>

class EvalLinks : public PHObject {

public:

  virtual ~EvalLinks() {}
  
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset() {return;}
  virtual int  isValid() const {return 0;}
  
  // container modifiers
  virtual void set_names(const std::string& left_name,
			 const std::string& right_name,
			 const std::string& weight_name) {return;}
  virtual void link(unsigned long long left_id, unsigned long long right_id, double weight) {return;}
  virtual void unlink(unsigned long long left_id, unsigned long long right_id) {return;}
  virtual void unlink_subleading() {return;}
  virtual void clear() {return;}

  // container status
  virtual size_t size() const {return 0;}
  virtual std::string get_name_left() const {return std::string();}
  virtual std::string get_name_right() const {return std::string();}
  virtual std::string get_name_weight() const {return std::string();} 
  virtual bool   has_link(unsigned long long left_id, unsigned long long right_id) const {return false;}
  virtual double get_weight(unsigned long long left_id, unsigned long long right_id) const {return NAN;}

  virtual bool   is_null_id(unsigned long long id) const {return true;}
  virtual unsigned long long get_null_id() const {return ULONG_LONG_MAX;}
  
  virtual std::set<unsigned long long> left(unsigned long long right_id) const {
    return std::set<unsigned long long>();
  }
  virtual std::set<unsigned long long> right(unsigned long long left_id) const {
    return std::set<unsigned long long>();
  }

  // access maximum weight association
  virtual unsigned long long max_left(unsigned long long right_id) const {return ULONG_LONG_MAX;}
  virtual unsigned long long max_right(unsigned long long left_id) const {return ULONG_LONG_MAX;}

protected:
  EvalLinks(const std::string& left_name,
	    const std::string& right_name,
	    const std::string& weight_name);
  
private:
  ClassDef(EvalLinks,1)
};

#endif // __EVALLINKS_H__
