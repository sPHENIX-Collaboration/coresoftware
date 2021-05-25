#ifdef __CINT__
//#pragma link off all globals;
//#pragma link off all classes;
//#pragma link off all functions;

// -! means no streamers which is root'ish for
// we don't want to save this class on a file
// (who wants to save a reconstruction module
// in a file - no point to that)
//#pragma link C++ class FieldSim-!;
#pragma link C++ class AnnularFieldSim;
#pragma link C++ class CylindricalFieldSim;
#pragma link C++ class MultiArray;

#endif /* __CINT__ */
