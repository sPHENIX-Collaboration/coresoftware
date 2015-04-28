#ifndef __BUFFER_CONSTANTS_H
#define __BUFFER_CONSTANTS_H

/* the header length values */ 
#define BUFFERHEADERLENGTH 4
#define EOBLENGTH 4
#define BUFFERSIZE (1024*1024)

#define BUFFERMARKER      0xffffffc0
#define GZBUFFERMARKER    0xfffffafe
#define LZO1XBUFFERMARKER 0xffffbbfe

#define BUFFERBLOCKSIZE 8192



#endif
