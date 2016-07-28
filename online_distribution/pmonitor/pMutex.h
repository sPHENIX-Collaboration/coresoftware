#include <pthread.h>

class pMutex { //++CINT
  
 public:
  
  pMutex(const int lockstatus = 0);
  virtual ~pMutex();

  virtual int Lock();
  virtual int tryLock();
  virtual int Release();

 private:

  pthread_mutex_t M;
  

};
