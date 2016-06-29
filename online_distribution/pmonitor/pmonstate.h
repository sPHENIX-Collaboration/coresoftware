#ifndef __PMONSTATE__
#define __PMONSTATE__

#define ETSTREAM   1001
#define FILESTREAM 1002
#define TESTSTREAM 1003
#define RCDAQSTREAM 1004

#include "Event/Eventiterator.h"

#include "Event/msg_profile.h"
#include "Event/msg_control.h"

#include <iostream>

class pmonstate
{ 

  friend std::ostream& operator<< (std::ostream& os, pmonstate &c)
    {
      //      std::cout << "Status:" << std::endl;

      if (c.isRunning())
        {
          std::cout << "Running at Event " << c.getNoevt() << "  ";
        }
      else
	{
        std::cout << "Not running  ";
	}

      if ( c.streamOpened() )
        {
          os << "Stream open:  ";
        }
      else
        {
          os << "no stream open.";
        }
      return os;
    };



 private:
  int stream;
  int initstatus;
  int eventloopexitstatus;
  int running;
  int number_of_events;
  char * Name;
  int idflag;

  Eventiterator * it;

 public:
  pmonstate()
    {
      eventloopexitstatus = 0;


      idflag = stream = running = number_of_events = 0;
      Name = 0;
      int pinit();
      msg_control *Message = new msg_control(MSG_TYPE_ONLINE,
                                         MSG_SOURCE_ET,
                                         MSG_SEV_INFORMATIONAL, "pmonitor");
      std::cout << *Message << "Welcome to pmonitor.  Type phelp() for help " << std::endl;
      delete Message;
      initstatus = pinit();
    };

  virtual ~pmonstate()
    {
      if (Name)
        delete [] Name;
    };

  inline virtual int streamOpened() const
    {
      return stream;
    };

  virtual int setETOpened(const char *name )
    {
      stream = ETSTREAM;
      Name = new char[strlen(name) + 1];
      strcpy (Name, name);
      return 0;
    };

  virtual int setFileOpened(const char *name)
    {
      stream = FILESTREAM;
      Name = new char[strlen(name) + 1];
      strcpy (Name, name);
      return 0;
    };

  virtual int setRCDAQOpened(const char *name)
    {
      stream = RCDAQSTREAM;
      Name = new char[strlen(name) + 1];
      strcpy (Name, name);
      return 0;
    };

  virtual int setTestOpened()
    {
      stream = TESTSTREAM;
      Name = new char[strlen("Test") + 1];
      strcpy (Name, "Test");
      return 0;
    };

  virtual int clearOpened()
    {
      stream = 0;
      delete [] Name;
      Name = 0;
      return 0;
    };


  virtual int isRunning() const
    {
      return running;
      return 0;
    };

  virtual int setRunning()
    {
      running = 1;
      return 0;
    };

  virtual int clearRunning()
    {
      running = 0;
      return 0;
    }

  virtual int isBroken() const
    {
      return initstatus;
    }


  virtual int getNoevt() const
    {
      return number_of_events;
    };

  virtual int setNoevt(const int n)
    {
      number_of_events = n ;
      return 0;
    };

  virtual int clearNoevt()
    {
      number_of_events = 0 ;
      return 0;
    };

  virtual int incrementNoevt()
    {
      number_of_events++;
      return 0;
    };

  virtual int getloopStatus()
    {
      return eventloopexitstatus;
    };

  virtual void setloopStatus(const int s)
    {
      eventloopexitstatus = s;
    };

  virtual int setIdentifyFlag(const int n = -1)
    {
      idflag = n;
      return 0;
    };

  virtual int isIdentifyFlag()
    {
      int i = idflag;
      if (idflag > 0)
        {
          idflag--;
        }
      return i;
    };

  virtual char *name() const
    {
      return Name;
    }
};

#endif /* __PMONSTATE__ */
