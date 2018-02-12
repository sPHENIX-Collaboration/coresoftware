#ifndef __TrkrDefUtil_H__
#define __TrkrDefUtil_H__


#ifdef __CINT__
#include <stdint.h>
#include <limits.h>
#else
#include <cstdint>
#include <climits>
#endif
#include <iostream>


/**
 * Define a namespace for typedefs
 */
namespace TrkrDefs
{
/// Key types
typedef uint32_t hitsetkey;  // 32 bit TrkrHitSet key type
typedef uint64_t cluskey;    // 64 but TrkrCluster id type
typedef uint32_t hitkey;     // 32 bit hit id type in TrkrCluster

// D. McGlinchey - I don't understand why this needs to be hidden from
//                 the dictionary generation ...
#ifndef __CINT__


/// Max values for keys (used as defaults or invalid values)
static hitsetkey HITSETKEYMAX = ULONG_MAX;
static cluskey CLUSKEYMAX = ULLONG_MAX;
static hitkey HITKEYMAX = ULONG_MAX;

#endif

/// Enumeration for tracker id to easily maintain consistency
enum TRKRID
{
  mvtx_id = 0,
  intt_id = 1,
  tpc_id = 2
};

}

/**
 * Utility class for parsing Trkr key types
 */
class TrkrDefUtil
{
public:

  /// ctor
  TrkrDefUtil() {};

  /// dtor
  ~TrkrDefUtil() {};

  /// Print the bits for each key type
  void print_bits(const TrkrDefs::hitsetkey key, std::ostream& os = std::cout);
  void print_bits(const TrkrDefs::cluskey key, std::ostream& os = std::cout);
  // void print_bits(const TrkrDefs::hitkey key, std::ostream& os = std::cout);

  /// Get the tracker ID from either key type
  uint8_t get_trackerid(const TrkrDefs::hitsetkey key);
  uint8_t get_trackerid(const TrkrDefs::cluskey key);

  /// Get the layer number from either key type
  uint8_t get_layer(const TrkrDefs::hitsetkey key);
  uint8_t get_layer(const TrkrDefs::cluskey key);

  /// Get the lower 32 bits for cluster keys only
  uint32_t get_index(const TrkrDefs::cluskey key);

  /// Get a valid low / hi range for hitsetkey given tracker id & layer
  TrkrDefs::hitsetkey get_hitsetkeylo(const TrkrDefs::TRKRID trkr_id);
  TrkrDefs::hitsetkey get_hitsetkeyhi(const TrkrDefs::TRKRID trkr_id);
  TrkrDefs::hitsetkey get_hitsetkeylo(const TrkrDefs::TRKRID trkr_id, const char lyr);
  TrkrDefs::hitsetkey get_hitsetkeyhi(const TrkrDefs::TRKRID trkr_id, const char lyr);

  /// Get a valid low / hi range for cluskey given tracker id & layer
  TrkrDefs::cluskey get_cluskeylo(const TrkrDefs::TRKRID trkr_id);
  TrkrDefs::cluskey get_cluskeyhi(const TrkrDefs::TRKRID trkr_id);
  TrkrDefs::cluskey get_cluskeylo(const TrkrDefs::TRKRID trkr_id, const char lyr);
  TrkrDefs::cluskey get_cluskeyhi(const TrkrDefs::TRKRID trkr_id, const char lyr);


protected:


  // hitsetkey layout:
  //  common upper 16 bits
  //   24 - 32  tracker id
  //   16 - 24  layer
  static const unsigned int bitshift_trackerid = 24; // 32 - 8
  static const unsigned int bitshift_layer = 16; // bitshift_trackerid - 8

  /// generate the common upper 16 bits for hitsetkey
  TrkrDefs::hitsetkey gen_hitsetkey(const TrkrDefs::TRKRID trkr_id, const char lyr);

  /// generate the common upper 16 bits for cluskey
  TrkrDefs::cluskey gen_cluskey(const TrkrDefs::TRKRID trkr_id, const char lyr);

};




#endif //__TrkrDefUtil_H__
