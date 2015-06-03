
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
  virtual void link(unsigned int left_id, unsigned int right_id, float weight) {return;}
  virtual void unlink(unsigned int left_id, unsigned int right_id) {return;}
  virtual void unlink_subleading() {return;}
  virtual void clear() {return;}

  // container status
  virtual size_t size() const {return 0;}
  virtual std::string get_name_left() const {return std::string();}
  virtual std::string get_name_right() const {return std::string();}
  virtual std::string get_name_weight() const {return std::string();} 
  virtual bool   has_link(unsigned int left_id, unsigned int right_id) const {return false;}
  virtual float  get_weight(unsigned int left_id, unsigned int right_id) const {return NAN;}

  virtual bool         is_null_id(unsigned int id) const {return true;}
  virtual unsigned int get_null_id() const {return UINT_MAX;}
  
  virtual std::set<unsigned int> left(unsigned int right_id) const {
    return std::set<unsigned int>();
  }
  virtual std::set<unsigned int> right(unsigned int left_id) const {
    return std::set<unsigned int>();
  }

  // access maximum weight association
  virtual unsigned int max_left(unsigned int right_id) const {return UINT_MAX;}
  virtual unsigned int max_right(unsigned int left_id) const {return UINT_MAX;}

protected:
  EvalLinks(const std::string& left_name,
	    const std::string& right_name,
	    const std::string& weight_name);
  
private:
  ClassDef(EvalLinks,1)
};

#endif // __EVALLINKS_H__
