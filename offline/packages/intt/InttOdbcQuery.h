#ifndef INTT_ODBC_QUERY_H
#define INTT_ODBC_QUERY_H

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>

namespace odbc
{
  class Statement;
}  // namespace odbc

class InttOdbcQuery
{
public:
  InttOdbcQuery() = default;
  ~InttOdbcQuery() = default;

  int Verbosity() {return m_verbosity;}
  int Verbosity(int verbosity) {return m_verbosity = verbosity;}

  int Query(int);

  bool IsStreaming() {return m_is_streaming;}
  std::vector<int> GetInttDACValues() {return m_intt_dac_values;} 
  const std::string &Type() {return m_type;}

private:
  int QueryStreaming(odbc::Statement *, int);
  int QueryType(odbc::Statement *, int);

  int QuerySingleDACValue(odbc::Statement *statement, int runnumber, int adc_value, int &DAC_value);
  int QueryAllDACValues(odbc::Statement *statement, int runnumber);

  static const int m_MAX_NUM_RETRIES = 3000;
  static const int m_MIN_SLEEP_DUR =  200; // milliseconds
  static const int m_MAX_SLEEP_DUR = 3000; // milliseconds

  int m_verbosity{0};
  bool m_query_successful{false};

  bool m_is_streaming{false};
  std::string m_type;
  std::array<std::set<std::string>, 8> m_file_set;

  std::vector<int> m_intt_dac_values;
  const std::vector<int> default_intt_dac_values{30, 45, 60, 90, 120, 150, 180, 210}; // note : most of physics runs have this setting.
};

#endif//INTT_ODBC_QUERY_H

