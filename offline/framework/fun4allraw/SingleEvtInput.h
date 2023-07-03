#ifndef FUN4ALLRAW_SINGLEEVTINPUT_H
#define FUN4ALLRAW_SINGLEEVTINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <list>
#include <map>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllEvtInputPoolManager;
class Packet;

class SingleEvtInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleEvtInput(const std::string &name, Fun4AllEvtInputPoolManager *inman);
  ~SingleEvtInput() override;
  Eventiterator *GetEventIterator() { return m_EventIterator; }
  void FillPool(const unsigned int nevents);
  int RunNumber() const { return m_RunNumber; }
  int fileopen(const std::string &filename) override;
  int fileclose() override;
  int AllDone() const {return m_AllDone;}
  void AllDone(const int i) {m_AllDone = i;}

 private:
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllEvtInputPoolManager *m_InputMgr = nullptr;
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  unsigned int m_EventNumberOffset = 1;  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber = 0;
  int m_EventsThisFile = 0;
  int m_AllDone = 0;
};

#endif
