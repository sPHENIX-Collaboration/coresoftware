#ifndef __PMONITOR__
#define __PMONITOR__

#include <Event/Event.h>

int pstatus ();                         

int ptestopen ();                       

int petopen (const char * etname);                       
int poncsetopen (const char * etname);                       

int poncsopen (const char * filename);  
int rcdaqopen (const char * ip=0);  
int pfileopen (const char * filename);  
int plistopen (const char * filename);  
int pstart (const int nevents);         
void prun ();                            
void prun (const int nevents);           
int pstart ();                          
int pstop ();                           
int pclose ();                          
//int pclear();
int pidentify (const int n);            
int pidentify ();                       
int pclearidentify ();                  
int pcontrol ();                        
int pcontrol (const int seconds);       
int pwait ();                           
int plock();                            
int prelease();                         
const char *pname();
int pexitstatus();                      
int pgui ();                            
int prmgui ();

int pinit ();
int process_event (Event * e);

int phsave (const char *filename ="histos.root"); 
void phelp();                          

#endif /* __PMONITOR__ */
