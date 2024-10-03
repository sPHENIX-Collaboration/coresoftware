#ifndef INTT_ODBC_QUERY_H
#define INTT_ODBC_QUERY_H

#include <array>
#include <set>
#include <string>

class InttOdbcQuery
{
public:
  InttOdbcQuery() = default;
  ~InttOdbcQuery() = default;

  int Verbosity() {return m_verbosity;}
  int Verbosity(int verbosity) {return m_verbosity = verbosity;}

  std::string Gl1Path() {return m_gl1_path;}
  std::string Gl1Path(std::string const& gl1_path) {return m_gl1_path = gl1_path;}

  std::string InttPath() {return m_intt_path;}
  std::string InttPath(std::string const& intt_path) {return m_intt_path = intt_path;}

  int Query(int);

  bool IsStreaming() {return m_is_streaming;}
  std::string Type() {return m_type;}

  std::set<std::string>::const_iterator Gl1FileListBegin() {return m_gl1_files.begin();}
  std::set<std::string>::const_iterator Gl1FileListEnd() {return m_gl1_files.end();}

  std::set<std::string>::const_iterator InttFileListBegin(int which_intt) {return m_intt_files.at(which_intt).begin();}
  std::set<std::string>::const_iterator InttFileListEnd(int which_intt) {return m_intt_files.at(which_intt).end();}

private:
  int QueryStreaming(void*, int);
  int QueryType(void*, int);
  int QueryGl1Files(void*, int);
  int QueryInttFiles(void*, int);

  int QueryFiles(void*, int, std::set<std::string>&, std::string const&, std::string const&); // helper method

  static const int m_MAX_NUM_RETRIES = 3000;
  static const int m_MIN_SLEEP_DUR =  200; // milliseconds
  static const int m_MAX_SLEEP_DUR = 3000; // milliseconds

  int m_verbosity{0};
  bool m_query_successful{false};
  std::string m_gl1_path =  "/sphenix/lustre01/sphnxpro/physics/GL1/";
  std::string m_intt_path = "/sphenix/lustre01/sphnxpro/physics/INTT/";

  bool m_is_streaming{false};
  std::string m_type;
  std::set<std::string> m_gl1_files;
  std::array<std::set<std::string>, 8> m_intt_files;
};

#endif//INTT_ODBC_QUERY_H

