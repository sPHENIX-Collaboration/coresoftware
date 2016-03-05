/*
**	Some standard definitions and typedefs for coding in the
**	PHENIX online system.
**
*/

#ifndef _phenixOnlineIncludeProtection_
#define _phenixOnlineIncludeProtection_



#include "phenixTypes.h"

#ifdef CONSTANT 
#undef CONSTANT
#endif
#ifdef __cplusplus
#define CONSTANT const
#else
#define CONSTANT static
#endif

#ifdef VXWORKS
#define INLINE_P inline
#define INLINE_D inline
#else
#define INLINE_P inline
#define INLINE_D inline
#endif

/*
**  Include malloc for some readback routines 
**  skip for DCM compiler.
*/
#if !defined(DCM) && !defined(VXWORKS)
#include <memory.h>
#endif

#include <assert.h>
#define Debug_Output

#define DWORD_SIZE sizeof(PHDWORD)
#define SWORD_SIZE sizeof(SWORD)

CONSTANT UINT DwordSize = DWORD_SIZE;
CONSTANT UINT SwordSize = SWORD_SIZE;

#define PRDF_BIG_ENDIAN 2
#define PRDF_LITTLE_ENDIAN 1

/*
**  We need to have all the various systems explicitly tested
**    with ifdef's to define the proper endianism.
*/
#define PRDF_LOCAL_ENDIANISM PRDF_LITTLE_ENDIAN

CONSTANT UINT PRDFbigEndian = PRDF_BIG_ENDIAN;
CONSTANT UINT PRDFlittleEndian = PRDF_LITTLE_ENDIAN;
CONSTANT UINT PRDFlocalEndianism = PRDF_LOCAL_ENDIANISM;

/*
** Function return values
*/
typedef int LOGIC_ret;
typedef UINT VALUE_ret;
typedef PHDWORD* PTR_ret;

#ifdef VXWORKS
#define TRUE  1
#define FALSE 0
#else
#ifndef TRUE
CONSTANT LOGIC_ret TRUE = 1;
#endif
#ifndef FALSE
CONSTANT LOGIC_ret FALSE = 0;
#endif
#endif

/*
**  Other constant definitions
*/
CONSTANT VALUE_ret valueFailure = 0xffffffff;
PTR_ret CONSTANT ptrFailure = 0;
CONSTANT LOGIC_ret logicFailure = FALSE;
CONSTANT LOGIC_ret logicSuccess = TRUE;
CONSTANT UINT maxByteValue = 0xff ;
CONSTANT UINT maxSwordValue = 0xffff ;
CONSTANT UINT maxDwordValue = 0xfffffffe ; /* reserve 0xffffffff for errors */


#define dwordCopy(out_ptr, in_ptr, numDwords) \
memcpy (out_ptr, in_ptr, 4*(numDwords))

#define dwordClear(out_ptr, numDwords) \
memset (out_ptr, 0, 4*(numDwords))

#define byteClear(out_ptr, numBytes) \
memset (out_ptr, 0, numBytes)

#define byteCopy(out_ptr, in_ptr, numBytes) \
memcpy (out_ptr, in_ptr, numBytes)

/*
**  Byte-swapping routines. There will likely be compiler or processor-specific
**    implementations of these.
*/
/*
**  A generic byte-swapping implementation that will likely work on any
**    machine but is likely to be slow. 
**
**  Note that the swap is not done "in-place".
*/
//inline int mlp();
//inline unsigned long singleDwordByteSwap(DWORD);
INLINE_D PHDWORD singleDwordByteSwap(PHDWORD inDword)
{
  PHDWORD outDword;

  outDword = (inDword & 0xFF) << 24;
  outDword |= (inDword >> 8 & 0xFF) << 16;
  outDword |= (inDword >> 16 & 0xFF) << 8;
  outDword |= (inDword >> 24 & 0xFF);
  return outDword;
}

/*
**  A routine to byte swap a sequence of dwords.
**
**  Since singleDwordByteSwap does not do an "in-place" swap, this routine
**  can be safely used with in_ptr = out_ptr.
*/
inline void dwordByteSwap (PHDWORD* out_ptr, PHDWORD* in_ptr, PHDWORD numDwords)
{
  UINT i; 
  for (i = 0; i< numDwords; i++) *out_ptr++ = singleDwordByteSwap(*in_ptr++);
}

/*
** Define macros to insert/extract bit fields from DWORD
*/

#define getWordMACRO(packet_ptr,offsetOfDWORD) (*((packet_ptr)+(offsetOfDWORD)))

#define getBitsMACRO(packet_ptr,offsetOfDWORD,offsetInDWORD,mask) \
             (((*((packet_ptr)+(offsetOfDWORD)))&(mask))>>offsetInDWORD)


#define setWordMACRO(packet_ptr,offsetOfDword,value) \
             (*(packet_ptr+offsetOfDword))=value


#define setBitsMACRO(packet_ptr,offsetOfDword,offsetInDword,mask,value) \
             (*(packet_ptr+offsetOfDword))&=(~mask); \
             (*(packet_ptr+offsetOfDword))|=(value<<offsetInDword)


#endif 

     /* end of ifndef _phenixOnlineIncludeProtection_ */
























