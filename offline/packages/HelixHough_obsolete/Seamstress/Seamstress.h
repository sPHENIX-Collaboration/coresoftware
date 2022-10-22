#ifndef SEAMSTRESS_SEAMSTRESS_H
#define SEAMSTRESS_SEAMSTRESS_H

#include <pthread.h>
#include <vector>

namespace SeamStress { class Needle; }

namespace SeamStress
{
  class Seamstress
  {
    public:
      Seamstress();
      Seamstress(const Seamstress &ss);
      ~Seamstress();
      
      static void *prepare(void *arg);
      
      void start();
      void stop();
      void sew();
      void rest();
      
      static void init_vector(unsigned long int N, std::vector<Seamstress> &vec);
      static std::vector<Seamstress*>* create_vector(unsigned long int N);
      
      
      Needle *needle;
      void *thread;
      
    private:
      bool gotime;
      bool end;
      bool queue_end;
      bool running;
      bool started;
      pthread_t pthread;
      pthread_attr_t attr;
      pthread_mutex_t mutex;
      pthread_mutexattr_t mattr;
      pthread_cond_t cond;
      pthread_cond_t waitcond;
  };
}


#endif
