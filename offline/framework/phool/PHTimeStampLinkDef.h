#ifdef __CINT__

#pragma link C++ class PHTimeStamp + ;

#ifndef __CLING__
#pragma link C++ function operator+ (const PHTimeStamp &, time_t);
#pragma link C++ function operator- (const PHTimeStamp &, time_t);
#pragma link C++ function operator- (const PHTimeStamp &, const PHTimeStamp &);
#pragma link C++ function operator<< (ostream &, const PHTimeStamp &);
#pragma link C++ function operator>> (istream &, const PHTimeStamp &);
#endif

#endif /* __CINT__ */
