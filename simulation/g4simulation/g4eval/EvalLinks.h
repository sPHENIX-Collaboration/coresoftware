
#ifndef __EVALLINKS_H__
#define __EVALLINKS_H__

#include <phool/PHObject.h>

#include <string>
#include <set>
#include <cmath>

class EvalLinks : public PHObject {

public:
  virtual ~EvalLinks() {}

  // container modifiers
  virtual void link(unsigned int left_id, unsigned int right_id, float purity) {return;}
  virtual void unlink(unsigned int left_id, unsigned int right_id) {return;}
  virtual void clear() {return;}

  // container status
  virtual void   print() const {return;}
  virtual size_t size() const {return 0;}
  virtual bool   has_link(unsigned int left_id, unsigned int right_id) const {return false;}
  virtual float  get_purity(unsigned int left_id, unsigned int right_id) const {return NAN;}

  virtual std::set<unsigned int> left(unsigned int right_id) const {
    return std::set<unsigned int>();
  }
  virtual std::set<unsigned int> right(unsigned int left_id) const {
    return std::set<unsigned int>();
  }

  // access maximum purity association
  virtual unsigned int max_left(unsigned int right_id) const {return 0xFFFFFFF;}
  virtual unsigned int max_right(unsigned int left_id) const {return 0xFFFFFFFF;}

protected:
  EvalLinks(const std::string left_name, std::string right_name) {}
  
private:
  ClassDef(EvalLinks,1)
};

#endif // __EVALLINKS_H__
