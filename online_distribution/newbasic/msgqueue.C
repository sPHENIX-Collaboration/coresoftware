
#include <msgqueue.h>

msgqueue::msgqueue (  char * token, const int id, 
	     int &status, const int inew )
{

  int ishmflg;

  tokenid = ftok ( token, id);
  wasnew = inew;
  if (inew)
    {
      ishmflg =  0666 | IPC_CREAT;
    }
  else
    {
      ishmflg =  0666;
    }
  msgid = msgget(tokenid,ishmflg);

  if (msgid < 0) 
    {
      status = msgid;
      return;
    }

  status = 0;
}
 
msgqueue::~msgqueue()
{
  if (wasnew)
    {

#if defined(SunOS) || defined(Linux)
      struct msqid_ds msgb;
      msgctl (msgid, IPC_RMID,&msgb);

#else
      msgctl (msgid, IPC_RMID);

#endif

    }
}

int msgqueue::receive( struct msgbuf *structure, int size)
{
  return msgrcv(msgid, structure, size,  0, 0);
}

int msgqueue::receive_nowait( struct msgbuf* structure, int size)
{
  return msgrcv(msgid, structure, size,  0, IPC_NOWAIT);
}

int msgqueue::send( struct msgbuf * structure, int size)
{
  return msgsnd(msgid, structure, size,  0);
}

int msgqueue::wait_for_ping()
{
  int i[2];
  return msgrcv(msgid, (struct msgbuf *) i, 8,  0, 0);
}

int msgqueue::send_ping()
{
  int i[2];
  i[0]= 1;
  i[1]= 0;

  return msgsnd(msgid, (struct msgbuf *) i , 8,  0);
}






