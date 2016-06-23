#include <pMutex.h>

pMutex::pMutex(const int lockstatus)
{
  pthread_mutex_init(&M, 0);

  if (lockstatus) pthread_mutex_lock(&M);

}

pMutex::~pMutex()
{
  pthread_mutex_destroy(&M);
}


  

int pMutex::Lock()
{
  return pthread_mutex_lock(&M);
}

int pMutex::tryLock()
{
  return pthread_mutex_trylock(&M);
}

int pMutex::Release()
{
  return pthread_mutex_unlock(&M);
}

