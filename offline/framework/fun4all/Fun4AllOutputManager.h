// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLOUTPUTMANAGER_H
#define FUN4ALL_FUN4ALLOUTPUTMANAGER_H

#include "Fun4AllBase.h"

#include <cstddef>  // for size_t
#include <string>
#include <vector>

class PHCompositeNode;

class Fun4AllOutputManager : public Fun4AllBase
{
 public:
  //! destructor
  ~Fun4AllOutputManager() override
  {
  }

  //! print method (dump event selector)
  void Print(const std::string &what = "ALL") const override;

  //! add a node in outputmanager
  virtual int AddNode(const std::string & /*nodename*/)
  {
    return 0;
  }

  //! add a runwise node in outputmanager
  virtual int AddRunNode(const std::string & /*nodename*/)
  {
    return 0;
  }

  //! not write a node in outputmanager
  virtual int StripNode(const std::string & /*nodename*/)
  {
    return 0;
  }

  //! not write a runwise node in outputmanager
  virtual int StripRunNode(const std::string & /*nodename*/)
  {
    return 0;
  }

  virtual void SaveRunNode(const int) { return; }

  /*! \brief
    add an event selector to the outputmanager.
    event will get written only if all event selectors process_event method
    return EVENT_OK
  */
  virtual int AddEventSelector(const std::string &recomodule);

  //! opens output file
  virtual int outfileopen(const std::string & /*nodename*/)
  {
    return 0;
  }

  //! Common method, called before calling virtual Write
  int WriteGeneric(PHCompositeNode *startNode);

  //! write starting from given node
  virtual int Write(PHCompositeNode * /*startNode*/)
  {
    return 0;
  }

  //! write specified node
  virtual int WriteNode(PHCompositeNode * /*thisNode*/)
  {
    return 0;
  }

  //! retrieves pointer to vector of event selector module names
  virtual std::vector<std::string> *EventSelector()
  {
    return &m_EventSelectorsVector;
  }

  //! retrieves pointer to vector of event selector module ids
  virtual std::vector<unsigned> *RecoModuleIndex()
  {
    return &m_RecoModuleIndexVector;
  }

  //! decides if event is to be written or not
  virtual int DoNotWriteEvent(std::vector<int> *retcodes) const;

  //! get number of Events
  virtual size_t EventsWritten() const { return m_NEvents; }
  //! increment number of events
  virtual void IncrementEvents(const unsigned int i) { m_NEvents += i; }
  //! get output file name
  virtual std::string OutFileName() const { return m_OutFileName; }
  void OutFileName(const std::string &name) { m_OutFileName = name; }

 protected:
  /*! 
    constructor.
    is protected since we do not want the  class to be created in root macros
  */
  Fun4AllOutputManager(const std::string &myname);
  Fun4AllOutputManager(const std::string &myname, const std::string &outfname);

 private:
  //! Number of Events
  unsigned int m_NEvents = 0;

  //! output file name
  std::string m_OutFileName;

  //! vector of event selectors modules
  std::vector<std::string> m_EventSelectorsVector;

  //! vector of associated module indexes
  std::vector<unsigned> m_RecoModuleIndexVector;
};

#endif
