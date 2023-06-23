#ifndef FUN4ALLRAW_SINGLEPRDFINPUT_H
#define FUN4ALLRAW_SINGLEPRDFINPUT_H

#include <fun4all/Fun4AllBase.h>

#include <list>
#include <map>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllPrdfInputPoolManager;
class Packet;

class SinglePrdfInput : public Fun4AllBase
{
 public:
  explicit SinglePrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman);
  ~SinglePrdfInput() override;
  Eventiterator *GetEventIterator() { return m_EventIterator; }
  void AddPrdfInputFile(const std::string &filename);
  void FillPool();
  int RunNumber() const { return m_RunNumber; }
  void UsedOneEvent() { m_PoolEvents--; }

 private:
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllPrdfInputPoolManager *m_InputMgr = nullptr;
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  unsigned int m_EventNumberOffset = 1;  // packet event counters start at 0 but we start with event number 1
  unsigned int m_PoolEvents = 0;
  unsigned int m_PoolDepth = 10;
  unsigned int m_LowWaterMark = 5;
  int m_RunNumber = 0;
  std::map<int, std::vector<Packet *>> m_EventMap;
  std::list<std::string> m_FileList;
  std::list<std::string> m_FileListCopy;
  std::list<std::string> m_FileListOpened;  // all files which were opened durin
};

#endif
