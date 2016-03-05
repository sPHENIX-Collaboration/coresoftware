       
#include "msg_control.h"

int  msg_control::xx_active =0;





msg_control::msg_control(const int mtype,
			 const int source,
			 const int severity,
			 const char * sourcecomponent)
{
  msg_type = mtype;
  msg_source = source;
  msg_severity = severity;
  storedseverity = severity;
  msg_sourcecomponent = new char[strlen(sourcecomponent) + 1 ];
  strcpy(msg_sourcecomponent,sourcecomponent);
}

msg_control::~msg_control(){
  delete [] msg_sourcecomponent;
};

int
msg_control::activate()
{
  int old = xx_active;
  xx_active = 1;
  return old;
}

int
msg_control::deactivate()
{
  int old = xx_active;
  xx_active = 0;
  return old;
}

void
msg_control::set_sourcecomponent (const char * msgsourcecomponent)
{

      delete msg_sourcecomponent;
      msg_sourcecomponent = new char[strlen(msgsourcecomponent) + 1 ];
      strcpy(msg_sourcecomponent,msgsourcecomponent);


}

void
msg_control::set_severity(const int severity) 
{
  storedseverity = msg_severity;
  msg_severity = severity; 
} // set_severity

OSTREAM&  operator<< (OSTREAM& os , msg_control &mc)
{
  if (msg_control::xx_active) {
    os << "<|A" 
       << SETW(MSG_TYPE_WIDTH) << mc.msg_type
       << " |B " 
       << SETW(MSG_SOURCE_WIDTH) << mc.msg_source
       << " |C " 
       << SETW(MSG_SEV_WIDTH) <<  mc.msg_severity
       << " |G "
       << SETW(strlen(mc.msg_sourcecomponent)+1) << mc.msg_sourcecomponent
       << " |>";
  }

  return os;
}

