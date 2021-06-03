#ifdef __CINT__
// -! means no streamers which is root'ish for
// we don't want to save this class on a file
// (who wants to save a reconstruction module
// in a file - no point to that)
//#pragma link C++ class FieldSim-!;
#pragma link C++ class Rossegger;

#endif /* __CINT__ */
