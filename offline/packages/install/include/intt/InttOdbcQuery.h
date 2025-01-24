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

  int Query(int);

  bool IsStreaming() {return m_is_streaming;}
  std::string Type() {return m_type;}
  std::set<std::string>::const_iterator FileListBegin(int which_intt) {return m_file_set.at(which_intt).begin();}
  std::set<std::string>::const_iterator FileListEnd(int which_intt) {return m_file_set.at(which_intt).end();}

private:
  int QueryStreaming(void*, int);
  int QueryType(void*, int);
  int QueryFiles(void*, int, int);

  static const int m_MAX_NUM_RETRIES = 3000;
  static const int m_MIN_SLEEP_DUR =  200; // milliseconds
  static const int m_MAX_SLEEP_DUR = 3000; // milliseconds

  int m_verbosity{0};
  bool m_query_successful{false};

  bool m_is_streaming{false};
  std::string m_type;
  std::array<std::set<std::string>, 8> m_file_set;
};

#endif//INTT_ODBC_QUERY_H

