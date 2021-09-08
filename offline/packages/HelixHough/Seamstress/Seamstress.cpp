#include "Seamstress.h"
#include "Needle.h"

#include <cstddef>
#include <memory>

using namespace std;


namespace SeamStress
{
  Seamstress::Seamstress()
  {
    gotime = false;
    end = false;
    queue_end = false;
    running = false;
    started = false;
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
    pthread_cond_init(&cond, NULL);
    pthread_cond_init(&waitcond, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  }


  Seamstress::Seamstress(const Seamstress &ss)
  {
    gotime = false;
    end = false;
    queue_end = false;
    running = false;
    started = false;
    pthread_mutexattr_init(&mattr);
    pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_NORMAL);
    pthread_mutex_init(&mutex, &mattr);
    pthread_cond_init(&cond, NULL);
    pthread_cond_init(&waitcond, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  }


  Seamstress::~Seamstress()
  {
    pthread_mutex_lock(&mutex);
    while(started==true)
    {
      pthread_cond_wait(&waitcond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mattr);
    pthread_cond_destroy(&cond);
    pthread_cond_destroy(&waitcond);
  }


  void *Seamstress::prepare(void *arg)
  {
    Seamstress *seamstress = (Seamstress*)arg;
    while((seamstress->end)==false)
    {
      pthread_mutex_lock(&(seamstress->mutex));
      while(seamstress->gotime==false)
      {
        pthread_cond_wait(&(seamstress->cond), &(seamstress->mutex));
      }
      if(seamstress->queue_end == true)
      {
        seamstress->end = true;
        pthread_mutex_unlock(&(seamstress->mutex));
        continue;
      }
      pthread_mutex_unlock(&(seamstress->mutex));
      
      
      seamstress->gotime = false;
      (*(seamstress->needle))(seamstress->thread);
      
      
      pthread_mutex_lock(&(seamstress->mutex));
      seamstress->running = false; 
      if(seamstress->queue_end==true){seamstress->end=true;}
      pthread_mutex_unlock(&(seamstress->mutex));
      pthread_cond_signal(&(seamstress->waitcond));
    }
    pthread_mutex_lock(&(seamstress->mutex));
    seamstress->started = false;
    seamstress->running = false;
    pthread_mutex_unlock(&(seamstress->mutex));
    pthread_cond_signal(&(seamstress->waitcond));
    pthread_exit(NULL);
  }


  void Seamstress::start()
  {
    pthread_mutex_lock(&mutex);
    if(started==false)
    {
      end = false;
      gotime = false;
      queue_end = false;
      running = false;
      started = true;
      pthread_mutex_unlock(&mutex);
      pthread_create(&pthread, &attr, &Seamstress::prepare, this);
    }
    else
    {
      pthread_mutex_unlock(&mutex);
    }
  }


  void Seamstress::stop()
  {
    pthread_mutex_lock(&mutex);
    queue_end = true;
    gotime = true;
    pthread_mutex_unlock(&mutex);
    pthread_cond_signal(&cond);
  }


  void Seamstress::sew()
  {
    pthread_mutex_lock(&mutex);
    gotime = true;
    running = true;
    pthread_mutex_unlock(&mutex);
    pthread_cond_signal(&cond);
  }


  void Seamstress::rest()
  {
    pthread_mutex_lock(&mutex);
    while(running==true)
    {
      pthread_cond_wait(&waitcond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
  }


  void Seamstress::init_vector(unsigned long int N, vector<Seamstress> &vec)
  {
    vec.clear();
    Seamstress st;
    
    for(unsigned int i=0;i<N;i++)
    {
      vec.push_back(st);
    }
    for(unsigned int i=0;i<N;i++)
    {
      vec[i].start();
    }
  }
  
  
  vector<Seamstress*>* Seamstress::create_vector(unsigned long int N)
  {
    vector<Seamstress*>* vec = new vector<Seamstress*>();
    for(unsigned long int i=0;i<N;i++){vec->push_back(new Seamstress());}
    for(unsigned long int i=0;i<N;i++){(*vec)[i]->start();}
    return vec;
  }
}


