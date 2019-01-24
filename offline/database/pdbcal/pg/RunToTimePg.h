#ifndef __RUNTOTIMEPG_HH__
#define __RUNTOTIMEPG_HH__

#include <pdbcalbase/RunToTime.h>
#include <map>
#include <string>

namespace odbc
{
class Connection;
}

class RunToTimePg : public RunToTime
{
 protected:
  RunToTimePg();

 public:
  virtual ~RunToTimePg();

  static RunToTimePg *instance()
  {
    return mySpecificCopy;
  }

  PHTimeStamp *getBeginTime(const int runNumber);
  PHTimeStamp *getEndTime(const int runNumber);
  int getRunNumber(const PHTimeStamp &ts);
  int DisconnectDB();
  static int Register();

 private:
  PHTimeStamp *getTime(const int runNumber, const std::string &what);
  int GetConnection();
  static RunToTimePg *mySpecificCopy;
  std::map<const int, PHTimeStamp *> beginruntimes;
  std::map<const int, PHTimeStamp *> endruntimes;
  odbc::Connection *con;
};

#endif /* __RUNTOTIMEPG_HH__ */
