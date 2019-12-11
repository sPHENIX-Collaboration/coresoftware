// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_RUNHEADERV1_H
#define FFAOBJECTS_RUNHEADERV1_H

#include "RunHeader.h"

#include <iostream>
#include <ctime>


class RunHeaderv1: public RunHeader
{
 public:
  RunHeaderv1();
  virtual ~RunHeaderv1() {}

  void Reset();
  void identify(std::ostream& os = std::cout) const;
  int isValid() const;

   int get_RunNumber() const {return RunNumber;}
   void set_RunNumber(const int run) {RunNumber= run; return;}

   double get_Bfield() const {return Bfield;}
   void set_Bfield(const double rval) {Bfield = rval; return;}

   time_t get_TimeStart() const {return 0;}
   void set_TimeStart(const time_t start);

   time_t get_TimeStop() const {return 0;}
   void set_TimeStop(const time_t stop);

 protected:
   int RunNumber;
   time_t TimeStart;
   time_t TimeStop;
   double Bfield;

   ClassDef(RunHeaderv1,1)

};

#endif /* __RUNHEADERV1_H */

