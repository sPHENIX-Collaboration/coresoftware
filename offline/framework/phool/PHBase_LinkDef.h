#ifdef __CINT__

#pragma link C++ typedef PHBoolean;
#pragma link C++ enum PHMessageType;
#pragma link C++ enum PHAccessType;
#pragma link C++ enum  PHTreeType;
#pragma link C++ function PHMessage(const std::string&, int, const std::string&);
#pragma link C++ class PHRandomSeed-!;
#pragma link C++ class PHFlag-! ;
#pragma link C++ class PHObject+ ;
#pragma link C++ class PHTimeServer-!;
#pragma link C++ class PHTimeStamp+;
#pragma link C++ class recoConsts-!;

#pragma link C++ function operator +  (const PHTimeStamp &, time_t);
#pragma link C++ function operator -  (const PHTimeStamp &, time_t);
#pragma link C++ function operator -  (const PHTimeStamp &, const PHTimeStamp &);
#pragma link C++ function operator << (ostream &, const PHTimeStamp &);
#pragma link C++ function operator >> (istream &, const PHTimeStamp &);

#endif /* __CINT__ */
