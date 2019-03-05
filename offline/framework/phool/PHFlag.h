#ifndef PHOOL_PHFLAG_H
#define PHOOL_PHFLAG_H

/*
  General purpose flag package:
  Flags are identified by their name, Print() prints them out sorted
  alphabetically.
  To create e.g. an int flag called MYFLAG with value 10 use
  set_IntFlag("MYFLAG",10);
  get_IntFlag("MYFLAG")  will return 10 now.
  If an unitialised flag is used you'll get a warning printed out,
  the return value in this case is 0 for Int, 0. for float and NULL for
  char *
*/

#include <map>
#include <string>

class PHFlag
{
 public:
  PHFlag() {}
  virtual ~PHFlag() {}
  virtual const std::string get_CharFlag(const std::string &name) const;
  virtual const std::string get_CharFlag(const std::string &name, const std::string &defaultval);
  virtual void set_CharFlag(const std::string &name, const std::string &flag);

  virtual double get_DoubleFlag(const std::string &name) const;
  virtual double get_DoubleFlag(const std::string &name, const double defaultval);
  virtual void set_DoubleFlag(const std::string &name, const double flag);

  virtual float get_FloatFlag(const std::string &name) const;
  virtual float get_FloatFlag(const std::string &name, const float defaultval);
  virtual void set_FloatFlag(const std::string &name, const float flag);

  virtual int get_IntFlag(const std::string &name) const;
  virtual int get_IntFlag(const std::string &name, const int defaultval);
  virtual void set_IntFlag(const std::string &name, const int flag);

  virtual void Print() const;
  virtual void PrintDoubleFlags() const;
  virtual void PrintIntFlags() const;
  virtual void PrintFloatFlags() const;
  virtual void PrintCharFlags() const;
  virtual void ReadFromFile(const std::string &name);
  virtual void WriteToFile(const std::string &name);

  virtual int FlagExist(const std::string &name) const;

  virtual const std::map<std::string, int> *IntMap() const { return &intflag; }
  virtual const std::map<std::string, float> *FloatMap() const { return &floatflag; }
  virtual const std::map<std::string, double> *DoubleMap() const { return &doubleflag; }
  virtual const std::map<std::string, std::string> *CharMap() const { return &charflag; }

 protected:
  std::map<std::string, int> intflag;
  std::map<std::string, double> doubleflag;
  std::map<std::string, float> floatflag;
  std::map<std::string, std::string> charflag;
};

#endif
