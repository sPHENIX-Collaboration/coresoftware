#ifndef __SUBEVT_CONSTANTS_H
#define __SUBEVT_CONSTANTS_H

/* the enum types for dump style */
#define EVT_DECIMAL     1
#define EVT_HEXADECIMAL 2
#define EVT_OCTAL       3

/* Misc. values  */
#define MAX_OUTLENGTH 80000

// the header length value
#define SEVTHEADERLENGTH 4



// ---------------------------------------------------------------------
//    IDCRAW requests the subevent to be copied without any decoding
#define IDCRAW   0 

// ---------------------------------------------------------------------
//    IDDGEN uses the standard decoding method imbedded in the subevent
//    header in the new data format
#define IDDGEN    1 

// ---------------------------------------------------------------------
//    IDHCPY requests only the subevent header (or the Event header) to be
//    copied:
#define IDHCPY    2 

// ---------------------------------------------------------------------
//    the next methods < 10 use what we consider standard methods by
//    now, i.e., no scheme proprietary to one particular hardware brand

#define ID1STR    3 
#define IDCSTR    4 
#define ID2EVT    5 
#define ID4EVT    6 
#define ID2SUP    7 

// ---------------------------------------------------------------------
// the next methods are for the hammond/g-2 board.

#define IDHAMMONDSET    31
#define IDHAMMOND       32

#define IDSAM           40

#define IDMIZNHC        41

#define IDDCFEM         51
#define IDTECFEM        52

#define IDSIS3300       55
#define IDCAENV792      56
#define IDCAENV785N     57

#define IDFIFOBOARD     58
#define IDRCPETDATA     59
#define IDBSPETDATA     60
#define IDUPPETDATA     61
#define IDUPPETDATA_V104     62
#define IDSIS3300R       65


#define IDSRSV01          70
#define IDUPPETPARAMS     71

#define IDDRS4V1        81
#define IDCAENV1742     85

#define IDPCONTAINER     89

#define IDFNALMWPC     90
#define IDFNALMWPCV2     91


#define IDDIGITIZERV1     92




#endif
