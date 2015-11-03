#ifndef __RUNTOTIME_HH__
#define __RUNTOTIME_HH__

class PHTimeStamp;

class RunToTime 
{

protected:

  RunToTime();
  virtual ~RunToTime();

  static  RunToTime *__instance; 
public:

  virtual PHTimeStamp *getBeginTime(const int runNumber) = 0;
  virtual PHTimeStamp *getEndTime(const int runNumber) = 0;
  virtual int getRunNumber(const PHTimeStamp& ts) = 0;
  virtual int DisconnectDB() = 0;
  static RunToTime *instance();
};


#endif /* __RUNTOTIME_HH__ */
