#ifndef SEAMSTRESS_PINCUSHION_H
#define SEAMSTRESS_PINCUSHION_H

#include "Seamstress.h"
#include "Needle.h"

namespace SeamStress
{
  template <class TClass>
  class Pincushion : public Needle
  {
    public:
      Pincushion<TClass>(TClass *glove, std::vector<Seamstress*> *ss)
      {
        hand = glove;
        seamstresses = ss;
      }
      virtual ~Pincushion<TClass>(){}
      
      
      void operator()(void *thread)
      {
        (hand->*pattern)(thread);
      }
      
      
      void sew(void (TClass::*design)(void*), std::vector<void*> threads)
      {
        pattern = design;
        unsigned long int ntot = threads.size();
        if(seamstresses->size() < ntot){ntot = seamstresses->size();}
        for(unsigned long int i=0;i<ntot;i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = threads[i];
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<ntot;i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewSoftly(void (TClass::*design)(void*))
      {
        pattern = design;
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewSoftly(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewStraight(void (TClass::*design)(void*))
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<seamstresses->size();i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewStraight(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<N;i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          (*seamstresses)[i]->sew();
        }
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void sewOpenlyStraight(void (TClass::*design)(void*))
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<seamstresses->size();i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlyStraight(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        std::vector<unsigned long int> quilt;
        for(unsigned long int i=0;i<N;i++){quilt.push_back(i);}
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = &(quilt[i]);
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlySoft(void (TClass::*design)(void*))
      {
        pattern = design;
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void sewOpenlySoft(void (TClass::*design)(void*), unsigned long int N)
      {
        pattern = design;
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->needle = this;
          (*seamstresses)[i]->thread = NULL;
          (*seamstresses)[i]->sew();
        }
      }
      
      
      void tieOff()
      {
        for(unsigned long int i=0;i<seamstresses->size();i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      void tieOff(unsigned long int N)
      {
        for(unsigned long int i=0;i<N;i++)
        {
          (*seamstresses)[i]->rest();
        }
      }
      
      
      std::vector<Seamstress*> *seamstresses;
      TClass *hand;
      void (TClass::*pattern)(void*);
  };
}


#endif
