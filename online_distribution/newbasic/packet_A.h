#ifndef __PACKET_A_H__
#define __PACKET_A_H__


#include "packet.h"
#include "packetPublic.h"
#include "packetConstants.h"
#include "packetRoutines.h"
#include "decoding_routines.h"

const char *get_mnemonic (const int structure,const int format);
const char *get_type_mnemonic (const int id);

/** in order to keep the top-level Packet class
    fully implementation-independent, we put 
    this Packet\_A class in, from which all current implementations
    derive.
*/

#ifndef __CINT__
class WINDOWSEXPORT Packet_A : public Packet
#else
class  Packet_A : public Packet
#endif
{

public:

  Packet_A(PACKET_ptr packet_ptr);
  Packet_A();
  ~Packet_A();

  // access to envelope information 
  int   getLength() const;
  //  int   getType() const;
  //  int   getDecoding() const;

  // some more header fields which are not yet implemented, marked "//*"
  //* int   gethdrVersion() const; // Version of header definition                    
  //* int   getHdrLength() const;     // inclusive of alignment data 
  //* int   getStatus() const;	       // Status bits describe frame errors/other 
  int   getErrorLength() const;    // Length of error block in Dwords
  int	getDebugLength() const;// Length of debug block in Dwords

  int	getIdentifier() const; // Identifier
  //* int	getEndianism() const;  // Big/little endian indicator
  int	getPadding() const;    // number of padding units

  int	getStructure() const;  // Structure of packet
  //* int	getWordSize() const;   // "Word" size used to store packet data
  //* int	getAddrLength() const; // number of bytes used for channel address
  //* int	getHitLength() const;  // Length of a single "hit" in bytes

  int	getHitFormat() const;  // Format of a single hit
  //* int	getNumEntries() const; // Number of "objects" stored in packet

  int	getDataLength() const;  // Format of a single hit


  // debugging-type information 


  void  identify( OSTREAM& =COUT) const;
  void  fullIdentify( OSTREAM& =COUT) const;

  void  dumpErrorBlock ( OSTREAM& =COUT ) ;
  void  dumpDebugBlock ( OSTREAM& =COUT ) ;


  // getting decoded values
  int    iValue(const int);
  int    iValue(const int,const char *);
  int    iValue(const int,const int);
  int    iValue(const int, const int, const char *){return 0;};
  int    iValue(const int channel, const int iy, const int iz) {return 0;};
  int    iValue(const int channel, const int iy, const int iz, const char *what) {return 0;};

  float  rValue(const int);
  float  rValue(const int,const char *);
  float  rValue(const int,const int);

  int    getArraylength(const char *);
  int    fillIntArray (int [], const int, int *,const char * what="");
  int    fillFloatArray (float [], const int, int *,const char * what="");
  int*   getIntArray (int *,const char * what="");
  float* getFloatArray (int *,const char *what="");

  // pointer or data based handling
  virtual int is_pointer_type() const;
  virtual int convert();

  int getCheckSumStatus() const;

  int copyMe(int dest [],  const int maxlength) const;

protected:

  //  PACKETHDR_ptr packetHdr;

  int    standardIntArray (int [], const int, int *,const char * what="");

  PACKET_ptr packet;  // storage for the packet


  int is_data_type;  // 0 is pointer based --  1 is data based

  int data1_length;
  int data2_length;
  int data3_length;
  int data4_length;
  int data5_length;
  int data6_length;
  int data7_length;
  int data8_length;
  int data9_length;
  int data10_length;

  int *decoded_data1;
  int *decoded_data2;
  int *decoded_data3;
  int *decoded_data4;
  int *decoded_data5;
  int *decoded_data6;
  int *decoded_data7;
  int *decoded_data8;
  int *decoded_data9;
  int *decoded_data10;

  virtual int *decode(int *) =0;

#ifdef LVL2_WINNT
  static void fix_endianess ( LONGLONG *x);
#else
  static void fix_endianess ( long long *x);
#endif

  static void fix_endianess ( double *x);
  static void fix_endianess ( char *str, const int length);

};



struct cdevIrData
{
      char m_irState[256];
      double m_tripletTrimCurrents[12];
      double m_irVacuum;
      double m_estBeamSizeYellowVert;
      double m_estBeamSizeYellowHorz;
      double m_estBeamSizeBlueVert;
      double m_estBeamSizeBlueHorz;
      double m_estimatedLuminosity;
      double m_betaStarYellowHorz;
      double m_betaStarBlueHorz;
      double m_betaStarYellowVert;
      double m_betaStarBlueVert;
      double m_polarPerBunchYellowX[360];
      double m_polarPerBunchYellowY[360];
      double m_polarPerBunchYellowZ[360];
      double m_polarPerBunchBlueX[360];
      double m_polarPerBunchBlueY[360];
      double m_polarPerBunchBlueZ[360];
      int m_avgOrbitDXBpmYellowHorzOdd;
      int m_avgOrbitDXBpmYellowHorzEven;
      int m_avgOrbitDXBpmYellowVertOdd;
      int m_avgOrbitDXBpmYellowVertEven;
      int m_avgOrbitDXBpmBlueHorzOdd;
      int m_avgOrbitDXBpmBlueHorzEven;
      int m_avgOrbitDXBpmBlueVertOdd;
      int m_avgOrbitDXBpmBlueVertEven;
      int m_experimentVertexX[100];
      int m_experimentVertexY[100];
      int m_experimentVertexZ[100];
      int m_vertexStartTime;
      int m_vertexEndTime;
      unsigned int m_datavalidMask; 

};

struct cdevRingData
{

  char   m_ringState[256]; // ejd91801
  char   m_ionSpecies[1024] ;
  double m_beamEnergy ;
  double m_gamma;
  int   m_stoneType;
  double m_momentumSpread;
  double m_bunchLength;
  int    m_intendedFillPattern[360];
  int    m_measuredFillPattern[360];
  double m_bunchOneRelativePhase;
  double m_synchrotronTune;
  double m_chromaticityVertical;
  double m_chromaticityHorizontal;
  int    m_polarizationFillPattern[360];
  int    m_timeOfFillStart;
  int    m_timeOfLuminosityStart;
  double m_emittanceVertical;
  double m_emittanceHorizontal;
  double m_betaIPMHorizontal;
  double m_betaIPMVertical;
  int    m_measuredPolarizationUp[360];
  int    m_measuredPolarizationDown[360];
  int    m_fillNumber;
  unsigned int m_datavalidMask;    // bit mask for data validity
  
};



struct cdevRingNoPolData
{

  char   m_ringState[256]; // ejd91801
  char   m_ionSpecies[1024] ;
  double m_beamEnergy ;
  double m_gamma;
  int   m_stoneType;
  double m_momentumSpread;
  double m_synchrotronTune;
  double m_chromaticityVertical;
  double m_chromaticityHorizontal;
  int    m_timeOfFillStart;
  int    m_timeOfLuminosityStart;
  double m_emittanceVertical;
  double m_emittanceHorizontal;
  double m_betaIPMHorizontal;
  double m_betaIPMVertical;
  int    m_fillNumber;
  unsigned int m_datavalidMask;    // bit mask for data validity
  
};

struct cdevBucketsData
{
  double m_bunchLength;
  double m_bunchOneRelativePhase;
  double m_fillPatternThreshold;
  int    m_intendedFillPattern[360];
  int    m_measuredFillPattern[360];
  int    m_polarizationFillPattern[360];
  unsigned int m_datavalidMask;    // bit mask for data validity
};

struct cdevRingPolData
{
  char   m_ringState[256]; // ejd91801
  char   m_ionSpecies[1024] ;
  double m_beamEnergy ;
  double m_gamma;
  int   m_stoneType;
  double m_momentumSpread;
  double m_bunchLength;
  int    m_intendedFillPattern[360];
  int    m_measuredFillPattern[360];
  double m_bunchOneRelativePhase;
  double m_synchrotronTune;
  double m_chromaticityVertical;
  double m_chromaticityHorizontal;
  int    m_polarizationFillPattern[360];
  int    m_timeOfFillStart;
  int    m_timeOfLuminosityStart;
  double m_emittanceVertical;
  double m_emittanceHorizontal;
  double m_betaIPMHorizontal;
  double m_betaIPMVertical;
  int    m_measuredPolarizationUp[360];
  int    m_measuredPolarization[360];
  unsigned int m_datavalidMask;    // bit mask for data validity
};
 
struct cdevWCMData
{
  int     cdevCaptureTimeStamp;
  //int dummy; //? 
  double  beamcurrent;
  float  bunchcurrent[360];

};

struct cdevDvmData
{
  double beamCurrent;
  double beamLifeTime;
};

struct cdevBPMData
{
  long    avgOrbTimeStamp;
  float  avgOrbPosition;
  float  avgOrbVariance;
  float  avgOrbStat;
  long   datavalidMask;
};

struct cdevMadchData
{
  int   cdevCaptureTimeStamp; //ejd 4/30/03 long to int
  double current;
};

struct cdevWCMHistory
{
  int counts;
  struct cdevWCMData reading[1];

};

struct cdevSISData
{

  int countRate;
  int totalCount;

};

struct cdevPolTargetData
{
  int positionEncLinear;
  int positionEncRot;
};

struct cdevPolarimeterData
{
  
  int   m_cdevCaptureTimeStamp;
  double runIdS;	// FILL.XXX --- where XXX is the run number 
  int startTimeS;	// Unix time
  int stopTimeS;	// Unix time
  char daqVersionS[80];
  char cutIdS[80];
  char targetIdS[80];	// "Horz.tagret3" or "Vert.target6" etc.
  int encoderPositionS[2];
  int statusS;	// bit pattern if <0 data is not usable
  char statusStringS[80];
  int totalCountsS;
  int upCountsS;
  int downCountsS;
  int unpolCountsS;
  int countsUpLeftS[360];
  int countsLeftS[360];
  int countsDownLeftS[360];
  int countsDownRightS[360];
  int countsRightS[360];
  int countsUpRightS[360];
  float avgAsymXS;
  float avgAsymX45S;
  float avgAsymX90S;    
  float avgAsymYS;
  float avgAsymErrorXS;
  float avgAsymErrorX45S;
  float avgAsymErrorX90S;
  float avgAsymErrorYS;
  float bunchAsymXS[360];
  float bunchAsymYS[360];
  float bunchAsymErrorXS[360];
  float bunchAsymErrorYS[360];
  float beamEnergyS;	// the same as ringSpec.color:beamEnergyM just for reference
  float analyzingPowerS;
  float analyzingPowerErrorS;
  int numberEventsS;	// provided by MCR before measurement
  int maxTimeS;	
  float polarizationM;
};



struct cdevPolarimeterZData
{
  
  int   m_cdevCaptureTimeStamp;
  double runIdS;	// FILL.XXX --- where XXX is the run number 
  int startTimeS;	// Unix time
  int stopTimeS;	// Unix time
  char daqVersionS[80];
  char cutIdS[80];
  char targetIdS[80];	// "Horz.tagret3" or "Vert.target6" etc.
  int encoderPositionS[2];
  int statusS;	// bit pattern if <0 data is not usable
  char statusStringS[80];
  int totalCountsS;
  int upCountsS;
  int downCountsS;
  int unpolCountsS;
  int countsUpLeftS[360];
  int countsLeftS[360];
  int countsDownLeftS[360];
  int countsDownRightS[360];
  int countsRightS[360];
  int countsUpRightS[360];
  float avgAsymXS;
  float avgAsymX45S;
  float avgAsymX90S;    
  float avgAsymYS;
  float avgAsymErrorXS;
  float avgAsymErrorX45S;
  float avgAsymErrorX90S;
  float avgAsymErrorYS;
  float bunchAsymXS[360];
  float bunchAsymYS[360];
  float bunchAsymErrorXS[360];
  float bunchAsymErrorYS[360];
  float beamEnergyS;	// the same as ringSpec.color:beamEnergyM just for reference
  float analyzingPowerS;
  float analyzingPowerErrorS;
  int numberEventsS;	// provided by MCR before measurement
  int maxTimeS;	
  float polarizationM;
};
#endif /* __PACKET_A_H__ */


