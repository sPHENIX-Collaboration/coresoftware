#ifndef __MSG_PROFILE_H__
#define __MSG_PROFILE_H__


// define some types of messages

#define MSG_TYPE_WIDTH       3
#define MSG_TYPE_MAX         7

#define MSG_TYPE_UNSPECIFIED 0
#define MSG_TYPE_ONLINE      1
#define MSG_TYPE_OFFLINE     2
#define MSG_TYPE_MONITORING  3
#define MSG_TYPE_CONTROL     4
#define MSG_TYPE_CODEDEBUG   5
#define MSG_TYPE_RUNTIME     6
#define MSG_TYPE_DEFAULT     MSG_TYPE_UNSPECIFIED

// define the source of a message

#define MSG_SOURCE_WIDTH        3
#define MSG_SOURCE_MAX          34

#define MSG_SOURCE_UNSPECIFIED  0
#define MSG_SOURCE_BEAMBEAM     1
#define MSG_SOURCE_BBC          1
#define MSG_SOURCE_MVD          2
#define MSG_SOURCE_DC           3
#define MSG_SOURCE_PC           4
#define MSG_SOURCE_TEC          5
#define MSG_SOURCE_RICH         6
#define MSG_SOURCE_TOF          7
#define MSG_SOURCE_PBSC         8
#define MSG_SOURCE_PBGL         9
#define MSG_SOURCE_MUTA        10
#define MSG_SOURCE_MUTC        11
#define MSG_SOURCE_MUID        12
#define MSG_SOURCE_HV          13         
#define MSG_SOURCE_ET          14
#define MSG_SOURCE_RC          15
#define MSG_SOURCE_EVB         16
#define MSG_SOURCE_ZDC         17
#define MSG_SOURCE_DAQMON      18
#define MSG_SOURCE_LVL1        19
#define MSG_SOURCE_LVL2        20
#define MSG_SOURCE_GL1         21
#define MSG_SOURCE_BUFFERBOX   22
#define MSG_SOURCE_AEROGEL     23
#define MSG_SOURCE_ERT         24
#define MSG_SOURCE_MPC         25
#define MSG_SOURCE_RXNP        26
#define MSG_SOURCE_LOCALPOL    27
#define MSG_SOURCE_MONITOR     28
#define MSG_SOURCE_MUTR        29
#define MSG_SOURCE_TOFW        30
#define MSG_SOURCE_CLOCK       31
#define MSG_SOURCE_VTX         32
#define MSG_SOURCE_FVTX        33

#define MSG_SOURCE_DEFAULT      MSG_SOURCE_UNSPECIFIED 

// define the severity of a message

#define MSG_SEV_WIDTH        2
#define MSG_SEV_MAX          5

#define MSG_SEV_INFORMATIONAL 0
#define MSG_SEV_WARNING      1
#define MSG_SEV_ERROR        2
#define MSG_SEV_SEVEREERROR  3
#define MSG_SEV_FATAL        4
#define MSG_SEV_DEFAULT      MSG_SEV_INFORMATIONAL

#endif 
