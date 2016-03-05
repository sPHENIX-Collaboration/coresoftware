#ifndef __PACKET_H__
#define __PACKET_H__


#include "generalDefinitions.h"

#include "event_io.h"

#define WINDOWSEXPORT

// --------------------------------------------------
// the virtual base base class for all Packets.
// --------------------------------------------------


/** 
This is the abstract packet class.
*/
#ifndef __CINT__
class WINDOWSEXPORT Packet
#else
class  Packet
#endif
{

  /* @name Description of the interface routines
     these are the interface routines to get at decoded values.
	iValue and rValue return int's and float's, respectively. 
	the optional "what" parameter is needed for some devices
	which have more than one type of information, such as the
	MIZAR board, which gets ADC and TDC values. 
	iValue(0), iValue(0,""), and iValue("ADC") all  
	return what we call "natural data", which means ADC information
	(I decided that), while iValue(0,"TDC") gives you the TDC
	value of channel 0. All Packet classes accept "RAW" as a type,
	directing them to return unprocessed/undecoded packet raw data.
  
	The array-type interface iValue(const int channel,const int iy)
	is just a convenient call for devices which inherently have
	two-dimensional arrays of data, such as the Hammond flash ADC
	board, where you get 8 channels with 512 time samples each.
	Another example is a group of several, say, LC2249W ADC's read out
	in one packet, so that iValue(0,2) refers to channel 0 of the
	second ADC in the group. 
  */




public:
  /// the virtual destructor
  inline virtual ~Packet() {};

  // **** getting decoded values ****



  /// iValue returns the value of a given channel as an int. 
  virtual int    iValue(const int channel) =0;

  /** with the "what" parameter you can decide which aspect 
      of the data you want to see (for devices which have more than one)
  */
  virtual int    iValue(const int channel, const char * what) =0;

  /** we have a few recent devices which have one more dimension 
      (such as card, time sample, channel)
  */
  virtual int    iValue(const int channel, const int y, const char * what) =0;

  /** this supports devices which are inherently organized as two-dimensional
      data, such as flash ADC's (channel vs time slice) 
  */
  virtual int    iValue(const int channel,const int iy) =0;

  /** this supports devices  organized as three-dimensional
      data (card vs channel vs time slice ) 
  */
  virtual int    iValue(const int channel,const int iy, const int iz) =0;

  /** this supports devices  organized as three-dimensional
      data (card vs channel vs time slice, with a "what" selection ) 
  */
  virtual int    iValue(const int channel,const int iy, const int iz, const char *what) =0;

  /** rValue returns the value of a given channel as a float */
  virtual float  rValue(const int channel) =0;

  /** dValue returns the value of a given channel as a double */
  virtual double  dValue(const int channel) 
    {return rValue(channel);};

  virtual double  dValue(const int channel, const char *what) 
    {return iValue(channel, what);};

  virtual double  dValue(const int channel, const int iy) 
    {return iValue(channel, iy);};

  /** lValue returns the value of a given channel as a long long */
  virtual long long  lValue(const int channel) 
    {return iValue(channel);};

  virtual long long  lValue(const int channel, const char *what) 
    {return iValue(channel,what);};

  virtual long long  lValue(const int channel, const int iy) 
    {return iValue(channel, iy);};


  /** with the "what" parameter you can decide which aspect 
      of the data you want to see (for devices which have more than one)
  */
  virtual float  rValue(const int channel, const char * what) =0;


  /** this supports devices which are inherently organized as two-dimensional
      data, such as flash ADC's (channel vs time slice) 
  */
  virtual float  rValue(const int channel, const int iy) =0;

  virtual void * pValue(const int chan)
    {
      return 0;
    }

  virtual void * pValue(const int chan, const char *what)
    {
      return 0;
    }

  virtual void * pValue(const int chan, const int iy)
    {
      return 0;
    }




  // *** now getting all the decoded values in one go ***
  // these routines get you all the decoded values into an array
  // of float's or int's. The what parameter has the usual meaning.

  /** getArraylength returns the length of the array needed to store the 
   decoded values.
  */
  virtual int    getArraylength(const char * what ="") =0;

  /** fillIntArray and fillFloatArray fill existing (user-supplied) arrays 
      with the decoded data
  */ 
  virtual int    fillIntArray (int destination[],    // the data go here 
			       const int length,      // space we have in destination
			       int * nw,              // words actually used
			       const char * what="") = 0; // type of data (see above)

  ///  fillFloatArray fills an array of floats
  virtual int    fillFloatArray (float destination[],    // the data go here 
			       const int length,      // space we have in destination
			       int * nw,              // words actually used
			       const char * what="") = 0; // type of data (see above)

  /** getIntArray and getFloatArray create a new array of the approriate size
      fill it with the decoded values, and return a pointer to the array. 
      nw is the length of the array created. 
  */
  virtual int*   getIntArray (int * nw,const char * ="") =0;

  ///  getFloatArray creates and returns an array of floats
  virtual float* getFloatArray (int * nw,const char * ="") =0;

  
  /// find out what type (pointer- or data based) packet object we have
  virtual int is_pointer_type() const = 0;

  /// convert from pointer- to data based object, if it is already data-based, do nothing.
  virtual int convert() =0;


  // access to envelope information:

  /** getLength() returns the length of the raw packet data. If you were to copy the data
      somewhere, the destination must be able to hold as many words.
  */
  virtual int   getLength() const = 0;

  //  virtual int   getType() const = 0;
  //  virtual int   getDecoding() const = 0;

  // some more header fields which are not yet implemented, marked "//*"
  //* virtual int   gethdrVersion() const = 0; // Version of header definition                    
  //* virtual int   getHdrLength() const = 0;     // inclusive of alignment data 
  //* virtual int   getStatus() const = 0;	       // Status bits describe frame errors/other 
  
  virtual int   getErrorLength() const = 0;    // Length of error block in Dwords

  /// get the length of the debug block
  virtual int	getDebugLength() const = 0;// Length of debug block in Dwords

  /// get the packet identifier
  virtual int	getIdentifier() const = 0; // Identifier

  //* virtual int	getEndianism() const = 0;  // Big/little endian indicator

  /// get the number of padding units in the packet data.
  virtual int	getPadding() const = 0;    // number of padding units

  /// get the structure of the packet data; unformatted, hitlist, etc.
  virtual int	getStructure() const = 0;  // Structure of packet

  //* virtual int	getWordSize() const = 0;   // "Word" size used to store packet data
  //* virtual int	getAddrLength() const = 0; // number of bytes used for channel address
  //* virtual int	getHitLength() const = 0;  // Length of a single "hit" in bytes

  /// get the hit format; in case of unformatted get the encoding scheme.
  virtual int	getHitFormat() const = 0;  // Format of a single hit
  //* virtual int	getNumEntries() const = 0; // Number of "objects" stored in packet

  /// get what the name says...
  virtual int	getDataLength() const = 0;  // Format of a single hit

  // debugging-type information
  // identify will write a short identification message to 
  // standard output or to the iOSTREAM provided. Nice to see.


  //#ifdef LVL2_WINNT
  ///see below for comments
  //  virtual void  identify( std::OSTREAM&os =std::COUT ) const = 0;
  //  virtual void  fullIdentify( std::OSTREAM& os =std::COUT ) const
  //    { identify(os);};
  //  virtual void  dump ( std::OSTREAM& =std::COUT )  = 0;
  //  virtual void  gdump (const int how = EVT_HEXADECIMAL, std::OSTREAM& =std::COUT) const = 0;


  /// write an identification message to the supplied OSTREAM.

  virtual void identify(std::ostream& os = std::cout) const = 0;



  /// write an indepth identification message to the supplied OSTREAM.
  virtual void fullIdentify(std::ostream& os = std::cout) const 
    { identify(os);};

  /**dump (either to standard output or the specified OSTREAM)
     the packet's data in some packet-specific way. 
  */

  virtual void dump(std::ostream& os = std::cout)  = 0;


  /** Since dump()
     requires the packet data to be consistent, gdump ("generic" dump)
     dumps the data without making assumptions about the integrity of
     the packet data. For many Packet classes dump just calls gdump
     because it is good enough. The "how" parameter for gdump specifies
     decimal, octal, or hexadecimal (0,1,2) dump.
  */

  virtual void gdump(const int how = EVT_HEXADECIMAL, std::ostream& os = std::cout) const  = 0;

  /** This functiion returns the status of the DCM transfer checksum.
      0  = cannot be calculated for this packet
      1  = ok
     -1  = fail

  */
  virtual int getCheckSumStatus() const { return 0;};
  virtual int   copyMe(int [],  const int maxlength) const { return 0;};


};

// some structures for the LVL triggers

// emcal 

struct emcChannelLongList
{
  int channel;
  int time;
  int highpost;
  int lowpost;
  int highpre;
  int lowpre;
};


#define RICH_TIME 0
#define RICH_POSTSAMPLE 1
#define RICH_PRESAMPLE 2

struct richChannelList
{
  int channel;
  int time;
  int post;
  int pre;
};

struct emcChannelShortList
{
  int channel;
  int gain;  // 0 = low, else high
  int time;
  int post;
  int pre;

};

struct tecChannelList
{
  int channel;
  int time;
  int value;

};


#endif /* __PACKET_H__ */
