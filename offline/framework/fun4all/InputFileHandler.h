#ifndef INPUTFILEHANDLER_H
#define INPUTFILEHANDLER_H

#include <list>
#include <string>

class InputFileHandler
{
 public:
  InputFileHandler() = default;
  virtual ~InputFileHandler() = default;
  virtual int fileopen(const std::string & /*filename*/) { return 0; }
  virtual int fileclose() { return -1; }
  int OpenNextFile();
  int AddListFile(const std::string &filename);
  int AddFile(const std::string &filename);
  void AddToFileOpened(const std::string &filename) { m_FileListOpened.push_back(filename); }
  void Print(const std::string &what = "ALL") const;
  int IsOpen() const { return m_IsOpen; }
  void IsOpen(const int i) { m_IsOpen = i; }
  void SetVerbosity(const int i) { m_Verbosity = i; }
  int GetVerbosity() const { return m_Verbosity; }
  void UpdateFileList();
  void FileName(const std::string &fn) { m_FileName = fn; }
  const std::string FileName() const { return m_FileName; }

 private:
  int m_IsOpen = 0;
  int m_Repeat = 0;
  int m_Verbosity = 0;
  std::string m_FileName;
  std::list<std::string> m_FileList;
  std::list<std::string> m_FileListCopy;
  std::list<std::string> m_FileListOpened;  // all files which were opened during running
};

#endif
