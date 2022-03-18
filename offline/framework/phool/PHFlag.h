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

  virtual double get_DoubleFlag(const std::string &name) const;
  virtual double get_DoubleFlag(const std::string &name, const double defaultval);
  virtual void set_DoubleFlag(const std::string &name, const double flag);

  virtual float get_FloatFlag(const std::string &name) const;
  virtual float get_FloatFlag(const std::string &name, const float defaultval);
  virtual void set_FloatFlag(const std::string &name, const float flag);

  virtual int get_IntFlag(const std::string &name) const;
  virtual int get_IntFlag(const std::string &name, const int defaultval);
  virtual void set_IntFlag(const std::string &name, const int flag);

  virtual uint64_t get_uint64Flag(const std::string &name) const;
  virtual uint64_t get_uint64Flag(const std::string &name, const uint64_t defaultval);
  virtual void set_uint64Flag(const std::string &name, const uint64_t flag);

  virtual const std::string get_StringFlag(const std::string &name) const;
  virtual const std::string get_StringFlag(const std::string &name, const std::string &defaultval);
  virtual void set_StringFlag(const std::string &name, const std::string &flag);

  virtual void Print() const;
  virtual void PrintDoubleFlags() const;
  virtual void PrintIntFlags() const;
  virtual void Printuint64Flags() const;
  virtual void PrintFloatFlags() const;
  virtual void PrintStringFlags() const;
  virtual void ReadFromFile(const std::string &name);
  virtual void WriteToFile(const std::string &name);

  virtual int FlagExist(const std::string &name) const;

  virtual const std::map<std::string, uint64_t> *uint64Map() const { return &m_UInt64FlagMap; }
  virtual const std::map<std::string, int> *IntMap() const { return &m_IntFlagMap; }
  virtual const std::map<std::string, float> *FloatMap() const { return &m_FloatFlagMap; }
  virtual const std::map<std::string, double> *DoubleMap() const { return &m_DoubleFlagMap; }
  virtual const std::map<std::string, std::string> *StringMap() const { return &m_StringFlagMap; }
  virtual void PrintStackTrace() const;

  virtual void ClearFlag(const std::string &name);
  virtual void ClearAll();

 protected:
  std::map<std::string, uint64_t> m_UInt64FlagMap;
  std::map<std::string, int> m_IntFlagMap;
  std::map<std::string, double> m_DoubleFlagMap;
  std::map<std::string, float> m_FloatFlagMap;
  std::map<std::string, std::string> m_StringFlagMap;
};

#endif
