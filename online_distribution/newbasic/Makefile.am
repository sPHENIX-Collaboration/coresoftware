AUTOMAKE_OPTIONS = foreign

if VXWORKS
include_HEADERS = oncsdefs.h \
  phenixTypes.h \
  EvtConstants.h \
  EvtStructures.h \
  oBuffer.h \
  EventTypes.h \
  BufferConstants.h 

else
# all definitions

allheaders = \
  ioselect.h \
  event_io.h \
  BoolDefs.h \
  devicedefs.h \
  oncsdefs.h \
  phenixOnline.h \
  phenixTypes.h \
  printmsg.h \
  Cframe.h \
  CframeV1.h \
  Cpacket.h \
  CpacketV1.h \
  dataBlock.h \
  dataBlockHdr.h \
  dataBlockPublic.h \
  errorBlock.h \
  fakeFrame.h \
  formatError.h \
  formatIO.h \
  frameHdr.h \
  frameHdrV1.h \
  framePackets.h \
  framePublic.h \
  frameRoutines.h \
  frameSourceIds.h \
  frameV1Public.h \
  packetHdr.h \
  packetHdrV1.h \
  packetPublic.h \
  packetRoutines.h \
  packetV1Public.h \
  A_Event.h \
  BufferConstants.h \
  Event.h \
  EventTypes.h \
  Eventiterator.h \
  EvtConstants.h \
  EvtStructures.h \
  simpleRandom.h \
  testEventiterator.h \
  fileEventiterator.h \
  listEventiterator.h \
  md5.h \
  PHmd5Utils.h \
  gen_utilities.h \
  dpipe_filter.h \
  PHmd5Value.h \
  buffer.h \
  prdfBuffer.h \
  decoding_routines.h \
  generalDefinitions.h \
  gzbuffer.h \
  lzobuffer.h \
  oBuffer.h \
  oEvent.h \
  ogzBuffer.h \
  olzoBuffer.h \
  oamlBuffer.h \
  oncsBuffer.h \
  oncsCollection.h \
  oncsEvent.h \
  oncsEventiterator.h \
  rcdaqEventiterator.h \
  oncsEvtConstants.h \
  oncsEvtStructures.h \
  oncsStructures.h \
  oncsSubConstants.h \
  oncsSub_id2evt.h \
  oncsSub_id4evt.h \
  oncsSub_idsis3300.h \
  oncsSub_idsis3300r.h \
  oncsSub_idcaenv792.h \
  oncsSub_idcaenv1742.h \
  oncsSub_idbspetdata.h \
  oncsSub_idrcpetdata.h \
  oncsSub_iduppetdata.h \
  oncsSub_iduppetdata_v104.h \
  oncsSub_iduppetparams.h \
  oncsSub_idsrs_v01.h \
  oncsSub_idfnalmwpc.h \
  oncsSub_idfnalmwpcv2.h \
  oncsSub_iddrs4v1.h \
  oncsSubevent.h \
  packet.h \
  packetConstants.h \
  packetHeaders.h \
  packet_A.h \
  packet_collection.h \
  packet_gl1.h \
  packet_gl1p.h \
  packet_gl1_evclocks.h \
  packet_gl1psum.h \
  packet_id2evt.h \
  packet_id4evt.h \
  packet_id4scaler.h \
  packet_w124.h \
  packet_cdevpolarimeter.h \
  packet_cdevir.h \
  packet_cdevwcm.h \
  packet_cdevbpm.h \
  packet_cdevmadch.h \
  packet_cdevpoltarget.h \
  packet_starscaler.h \
  packet_hbd_fpga.h \
  packet_hbd_fpgashort.h \
  packet_cdevring.h \
  packet_cdevringnopol.h \
  packet_cdevbuckets.h \
  packet_cdevdvm.h \
  packet_idcdevdescr.h


allsources = \
  Cframe.C \
  CframeV1.C \
  Cpacket.C \
  CpacketV1.C \
  dataBlock.C \
  errorBlock.C \
  fakeFrame.C \
  formatError.C \
  formatIO.C \
  framePackets.C \
  frameRoutines.C \
  packetRoutines.C \
  A_Event.cc \
  buffer.cc \
  prdfBuffer.cc \
  decoding_routines.cc \
  evt_mnemonic.cc \
  simpleRandom.cc \
  testEventiterator.cc \
  fileEventiterator.cc \
  listEventiterator.cc \
  md5.cc \
  PHmd5Utils.cc \
  PHmd5Value.cc \
  gzbuffer.cc \
  lzobuffer.cc \
  oBuffer.cc \
  oEvent.cc \
  ogzBuffer.cc \
  olzoBuffer.cc \
  oamlBuffer.cc \
  oncsBuffer.cc \
  oncsEvent.cc \
  oncsEventiterator.cc \
  rcdaqEventiterator.cc \
  oncsSub_idcstr.cc \
  oncsSub_id2evt.cc \
  oncsSub_id4evt.cc \
  oncsSub_idsis3300.cc \
  oncsSub_idsis3300r.cc \
  oncsSub_idcaenv792.cc \
  oncsSub_idcaenv1742.cc \
  oncsSub_idbspetdata.cc \
  oncsSub_idrcpetdata.cc \
  oncsSub_iduppetdata.cc \
  oncsSub_iduppetdata_v104.cc \
  oncsSub_iduppetparams.cc \
  oncsSub_idsrs_v01.cc \
  oncsSub_idfnalmwpc.cc \
  oncsSub_idfnalmwpcv2.cc \
  oncsSub_iddrs4v1.cc \
  oncsSubevent.cc \
  oncs_mnemonic.cc \
  packet_A.cc \
  packet_gl1.cc \
  packet_gl1p.cc \
  packet_gl1_evclocks.cc \
  packet_gl1psum.cc \
  packet_idcstr.cc \
  packet_id2evt.cc \
  packet_id4evt.cc \
  packet_id4scaler.cc \
  packet_mnemonic.cc \
  packet_w124.cc \
  packet_cdevpolarimeter.cc \
  packet_cdevir.cc \
  packet_cdevwcm.cc \
  packet_cdevbpm.cc \
  packet_cdevmadch.cc \
  packet_cdevpoltarget.cc \
  packet_starscaler.cc \
  packet_hbd_fpga.cc \
  packet_hbd_fpgashort.cc \
  packet_cdevring.cc \
  packet_cdevringnopol.cc \
  packet_cdevbuckets.cc \
  packet_cdevdvm.cc \
  packet_idcdevdescr.cc 

AM_LDFLAGS = -L$(libdir) -L$(ONLINE_MAIN)/lib -L$(OPT_SPHENIX)/lib

BUILT_SOURCES = ioselect.h 


msgheaders = \
  msg_control.h \
  msg_profile.h \
  msg_buffer.h \
  date_filter_msg_buffer.h \
  filter_msg_buffer.h \
  remote_msg_buffer.h \
  strnstr.h 

INCLUDES = -I$(includedir)  -I@ROOTINC@ -I$(OPT_SPHENIX)/include

noinst_HEADERS = \
  eventLinkDef.h \
  all.h

include_HEADERS = $(allheaders) $(msgheaders)

CLEANFILES = \
  ioselect.h \
  event_dict.C \
  event_dict.h \
  setup.com \
  msg_dict.C \
  msg_dict.h

lib_LTLIBRARIES = \
  libmessage.la \
  libRootmessage.la \
  libNoRootEvent.la \
  libEvent.la


bin_PROGRAMS = \
  dpipe \
  ddump \
  dlist \
  changeid \
  changehitformat \
  eventcombiner \
  prdf2prdf \
  prdfcheck \
  prdfsplit 


dpipe_SOURCES = dpipe.cc
dlist_SOURCES = dlist.cc
ddump_SOURCES = ddump.cc
eventcombiner_SOURCES = eventcombiner.cc 

if LINUX
dpipe_LDFLAGS = -Wl,-export-dynamic
endif

changeid_SOURCES = changeid.cc
changehitformat_SOURCES	= changehitformat.cc
prdf2prdf_SOURCES = prdf2prdf.cc
prdfcheck_SOURCES = prdfcheck.cc
prdfsplit_SOURCES = prdfsplit.cc


dpipe_LDADD = libNoRootEvent.la libmessage.la  -ldl  @LZOLIB@
dlist_LDADD = libNoRootEvent.la libmessage.la  -ldl  @LZOLIB@
ddump_LDADD = libNoRootEvent.la libmessage.la  -ldl @LZOLIB@
eventcombiner_LDADD = libNoRootEvent.la libmessage.la -ldl  @LZOLIB@


changeid_LDADD = libNoRootEvent.la libmessage.la -ldl @LZOLIB@
changehitformat_LDADD = libNoRootEvent.la libmessage.la -ldl @LZOLIB@
prdf2prdf_LDADD = libNoRootEvent.la libmessage.la -ldl  @LZOLIB@
prdfcheck_LDADD = libNoRootEvent.la -ldl  @LZOLIB@
prdfsplit_LDADD = libNoRootEvent.la -ldl  @LZOLIB@


libmessage_la_SOURCES = \
  date_filter_msg_buffer.cc \
  filter_msg_buffer.cc \
  remote_msg_buffer.cc \
  msg_buffer.cc \
  msg_control.cc \
  strnstr.cc 


libEvent_la_SOURCES =  event_dict.C 
libEvent_la_LIBADD = libNoRootEvent.la libRootmessage.la  @ROOTGLIBS@  -lz @LZOLIB@

libNoRootEvent_la_SOURCES = $(allsources) 
libNoRootEvent_la_LIBADD = libmessage.la  -lz @LZOLIB@


# because this if statement contains dependencies, no more definitions after
# this

if NEWIO
ioselect.h:
	echo "// *** this is a generated file. Do not commit, do not edit" > ioselect.h
	echo "#ifndef __IOSELECT_H__" >> ioselect.h
	echo "#define __IOSELECT_H__" >> ioselect.h
	echo  >> ioselect.h
	echo "#define STREAMBUF_NEW_IOSTREAM" >> ioselect.h
	echo  >> ioselect.h
	echo "#endif" >> ioselect.h
else
ioselect.h:
	echo "// *** this is a generated file. Do not commit, do not edit" > ioselect.h
	echo "#ifndef __IOSELECT_H__" >> ioselect.h
	echo "#define __IOSELECT_H__" >> ioselect.h
	echo  >> ioselect.h
	echo "#endif" >> ioselect.h

endif

if LINUX

libRootmessage_la_SOURCES = \
  msg_dict.C

msg_dict.C : \
  msg_control.h \
  date_filter_msg_buffer.h \
  filter_msg_buffer.h \
  remote_msg_buffer.h \
  msgLinkDef.h
	rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^

libRootmessage_la_LIBADD = libmessage.la

else

libRootmessage_la_SOURCES = \
  date_filter_msg_buffer.cc \
  filter_msg_buffer.cc \
  msg_buffer.cc \
  msg_control.cc \
  remote_msg_buffer.cc \
  strnstr.cc

endif



event_dict.C : \
  Event.h  \
  Eventiterator.h \
  fileEventiterator.h \
  listEventiterator.h \
  oncsEventiterator.h \
  rcdaqEventiterator.h \
  packet.h \
  testEventiterator.h \
  eventLinkDef.h
	rootcint -f $@ -c $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $^



# clean cache dir on Solaris 5.8
clean-local:
	rm -rf SunWS_cache

#if VXWORKS
endif
