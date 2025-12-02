// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLHISTOMANAGER_H
#define FUN4ALL_FUN4ALLHISTOMANAGER_H

#include "Fun4AllBase.h"

#include <map>
#include <string>

class Fun4AllOutputManager;
class TNamed;

class Fun4AllHistoManager : public Fun4AllBase
{
 public:
  explicit Fun4AllHistoManager(const std::string &name);
  ~Fun4AllHistoManager() override;

  void Print(const std::string &what = "ALL") const override;

  //! Register histogram or TTree object
  //! For histograms, enforce error calculation and propagation
  bool registerHisto(const std::string &hname, TNamed *h1d, const int replace = 0);

  //! Register histogram or TTree object
  //! For histograms, enforce error calculation and propagation
  bool registerHisto(TNamed *h1d, const int replace = 0);

  template <typename T>
  T *makeHisto(T *t)
  {
    if (not registerHisto(t))
    {
      delete t;
      t = nullptr;
    }
    return t;
  }
  int isHistoRegistered(const std::string &name) const;
  TNamed *getHisto(const std::string &hname) const;
  TNamed *getHisto(const unsigned int ihisto) const;
  std::string getHistoName(const unsigned int ihisto) const;
  unsigned int nHistos() const { return Histo.size(); }
  void Reset();
  int RunAfterClosing();
  int dumpHistos(const std::string &filename = "", const std::string &openmode = "RECREATE");
  const std::string &OutFileName() { return m_outfilename; }
  void setOutfileName(const std::string &filename) { m_outfilename = filename; }
  void SetClosingScript(const std::string &script) { m_RunAfterClosingScript = script; }
  void SetClosingScriptArgs(const std::string &args) { m_ClosingArgs = args; }
  void segment(const int segment) { m_CurrentSegment = segment; }
  int GetEventNumberRollover() const { return m_EventRollover; }
  void SetEventNumberRollover(const int evtno) { m_EventRollover = evtno; }
  int LastEventNumber() const { return m_LastEventNumber; }
  void SetLastEventNumber(int ival) { m_LastEventNumber = ival; }
  void UpdateLastEvent() { m_LastEventNumber += m_EventRollover; }
  void InitializeLastEvent(int eventnumber);
  void StartSegment(int iseg) { m_CurrentSegment = iseg; }
  void UseFileRule(bool b = true) { m_UseFileRuleFlag = b; }
  bool ApplyFileRule() const { return m_UseFileRuleFlag; }
  void CopyRolloverSetting(const Fun4AllOutputManager *outman);
  const std::string &LastClosedFileName() const { return m_LastClosedFileName; }
  bool isEmpty() const;

private:
  bool m_LastEventInitializedFlag{false};
  bool m_UseFileRuleFlag{false};
  int m_CurrentSegment{0};
  int m_EventRollover{0};
  int m_LastEventNumber{std::numeric_limits<int>::max()};

  std::string m_FileRule{"-%08d-%05d"};
  std::string m_outfilename;
  std::string m_RunAfterClosingScript;
  std::string m_ClosingArgs;
  std::string m_LastClosedFileName;
  std::map<const std::string, TNamed *> Histo;
};

#endif /* __FUN4ALLHISTOMANAGER_H */
