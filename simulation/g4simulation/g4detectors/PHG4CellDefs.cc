#include "PHG4CellDefs.h"

#include <phool/phool.h>

#include <iostream>

unsigned short
generic_lower_16bit_key(const PHG4CellDefs::keytype key, const PHG4CellDefs::CellBinning binning);

unsigned short
generic_upper_16bit_key(const PHG4CellDefs::keytype key, const PHG4CellDefs::CellBinning binning);

PHG4CellDefs::keytype
generic_16bit_genkey(const unsigned short detid, const PHG4CellDefs::CellBinning binning, const unsigned short upper16bits, const unsigned short lower16bits);

PHG4CellDefs::keytype
generic_32bit_genkey(const unsigned short detid, const PHG4CellDefs::CellBinning binning, const unsigned int bit32);

using namespace std;

PHG4CellDefs::keytype
PHG4CellDefs::SizeBinning::genkey(const unsigned short detid, const unsigned short zbin, const unsigned short iphi)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(detid, sizebinning, zbin, iphi);
  return key;
}

unsigned short int
PHG4CellDefs::SizeBinning::get_phibin(const PHG4CellDefs::keytype key)
{
  unsigned short int phibin = generic_lower_16bit_key(key, sizebinning);
  return phibin;
}

unsigned short int
PHG4CellDefs::SizeBinning::get_zbin(const PHG4CellDefs::keytype key)
{
  unsigned short int zbin = generic_upper_16bit_key(key, sizebinning);
  return zbin;
}

PHG4CellDefs::keytype
PHG4CellDefs::EtaPhiBinning::genkey(const unsigned short detid, const unsigned short iphi, const unsigned short ieta)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(detid, etaphibinning, iphi, ieta);
  return key;
}

unsigned short int
PHG4CellDefs::EtaPhiBinning::get_etabin(const PHG4CellDefs::keytype key)
{
  unsigned short int etabin = generic_lower_16bit_key(key, etaphibinning);
  return etabin;
}

unsigned short int
PHG4CellDefs::EtaPhiBinning::get_phibin(const PHG4CellDefs::keytype key)
{
  unsigned short int phibin = generic_upper_16bit_key(key, etaphibinning);
  return phibin;
}

PHG4CellDefs::keytype
PHG4CellDefs::SpacalBinning::genkey(const unsigned short etabin, const unsigned short phibin, const unsigned short fiberid)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(etabin, spacalbinning, phibin, fiberid);
  return key;
}

unsigned short int
PHG4CellDefs::SpacalBinning::get_etabin(const PHG4CellDefs::keytype key)
{
  unsigned long long tmp = key >> 48;
  unsigned short int etabin = tmp;
  return etabin;
}

unsigned short int
PHG4CellDefs::SpacalBinning::get_phibin(const PHG4CellDefs::keytype key)
{
  unsigned short int phibin = generic_upper_16bit_key(key, spacalbinning);
  return phibin;
}

unsigned short int
PHG4CellDefs::SpacalBinning::get_fiberid(const PHG4CellDefs::keytype key)
{
  unsigned short int fiberid = generic_lower_16bit_key(key, spacalbinning);
  return fiberid;
}


PHG4CellDefs::keytype
PHG4CellDefs::ScintillatorSlatBinning::genkey(const unsigned short detid, const unsigned short icolumn, const unsigned short irow)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(detid, scintillatorslatbinning, icolumn, irow);
  return key;
}

unsigned short int
PHG4CellDefs::ScintillatorSlatBinning::get_row(PHG4CellDefs::keytype key)
{
  unsigned short int rowbin = generic_lower_16bit_key(key, scintillatorslatbinning);
  return rowbin;
}

unsigned short int
PHG4CellDefs::ScintillatorSlatBinning::get_column(PHG4CellDefs::keytype key)
{
  unsigned short int columnbin = generic_upper_16bit_key(key, scintillatorslatbinning);
  return columnbin;
}
PHG4CellDefs::keytype
PHG4CellDefs::EtaXsizeBinning::genkey(const unsigned short detid, const unsigned short ixbin, const unsigned short ieta)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(detid, etaxsizebinning, ixbin, ieta);
  return key;
}

unsigned short int
PHG4CellDefs::EtaXsizeBinning::get_etabin(const PHG4CellDefs::keytype key)
{
  unsigned short int etabin = generic_lower_16bit_key(key, etaxsizebinning);
  return etabin;
}

unsigned short int
PHG4CellDefs::EtaXsizeBinning::get_xsizebin(const PHG4CellDefs::keytype key)
{
  unsigned short int etabin = generic_upper_16bit_key(key, etaxsizebinning);
  return etabin;
}

PHG4CellDefs::keytype
PHG4CellDefs::MapsBinning::genkey(const unsigned short detid, const unsigned int bit32_index)
{
  PHG4CellDefs::keytype key = generic_32bit_genkey(detid, mapsbinning, bit32_index);
  return key;
}

PHG4CellDefs::keytype
PHG4CellDefs::TPCBinning::genkey(const unsigned short detid, const unsigned short mod, const unsigned short pad)
{
  PHG4CellDefs::keytype key = generic_16bit_genkey(detid, tpcbinning, mod, pad);
  return key;
}

bool
PHG4CellDefs::has_binning(PHG4CellDefs::keytype key, PHG4CellDefs::CellBinning binning)
{
 keytype tmp = (key >> bitshift_binning) & 0xFFFF;
 if (tmp == binning)
   {
     return true;
   }
 return false;
}

short
PHG4CellDefs::get_binning(PHG4CellDefs::keytype key)
{
  keytype tmp = (key >> bitshift_binning) & 0xFFFF;
  short int i = tmp;
  return i;
}

short int
PHG4CellDefs::get_detid(const PHG4CellDefs::keytype key)
{
  keytype tmp = (key >> bitshift_layer);
  return tmp;
}

unsigned short 
generic_lower_16bit_key(const PHG4CellDefs::keytype key, const PHG4CellDefs::CellBinning binning)
{
  // check correct binning first
 PHG4CellDefs::keytype tmp = binning;
 tmp = (tmp << PHG4CellDefs::bitshift_binning);
 if ((key & tmp) == tmp)
   {
     unsigned short int low16bitkey = (key & 0xFFFF);
     return low16bitkey;
   }
 cout << PHWHERE << " could not decode 0x" << hex <<  key << dec << endl; 
 cout << "key 0x" << hex << key << ", binning: 0x" << tmp
      << " and: " <<  (key & tmp) << dec << endl;
 exit(1);
}

unsigned short 
generic_upper_16bit_key(const PHG4CellDefs::keytype key, const PHG4CellDefs::CellBinning binning)
{
  // check correct binning first
 PHG4CellDefs::keytype tmp = binning;
 tmp = (tmp << PHG4CellDefs::bitshift_binning);
 if ((key & tmp) == tmp)
   {
     PHG4CellDefs::keytype keytmp = key >> PHG4CellDefs::bitshift_upperkey;
     unsigned short int hi16bitkey = (keytmp & 0xFFFF);
     return hi16bitkey;
   }
 cout << PHWHERE << " could not decode 0x" << hex <<  key << dec << endl; 
 exit(1);
}

PHG4CellDefs::keytype
generic_16bit_genkey(const unsigned short detid, const PHG4CellDefs::CellBinning binning, const unsigned short upper16bits, const unsigned short lower16bits)
{
  PHG4CellDefs::keytype tmp = detid;
  PHG4CellDefs::keytype key = tmp << PHG4CellDefs::bitshift_layer; // layer/detector id used by extrating ranges
  tmp = binning;
  key |= (tmp << PHG4CellDefs::bitshift_binning); // binning method used to decode the key
  tmp = upper16bits;
  key |= (tmp << PHG4CellDefs::bitshift_upperkey); // upper bits used by column, so we can easily extract 
                                  // slats by column which are combined to towers
  key |= lower16bits;
  return key;
}

PHG4CellDefs::keytype
generic_32bit_genkey(const unsigned short detid, const PHG4CellDefs::CellBinning binning, const unsigned int bit32)
{
  PHG4CellDefs::keytype tmp = detid;
  PHG4CellDefs::keytype key = tmp << PHG4CellDefs::bitshift_layer; // layer/detector id used by extrating ranges
  tmp = binning;
  key |= (tmp << PHG4CellDefs::bitshift_binning); // binning method used to decode the key
  key |= bit32;
  return key;
}

