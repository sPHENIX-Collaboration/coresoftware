#ifndef __SUBEVTSTRUCTURES_H__
#define __SUBEVTSTRUCTURES_H__


typedef struct subevt_data
 {
   int   sub_length;
   short sub_id;
   short sub_type;
   short sub_decoding;
   short sub_padding;
   short reserved[2];
   int   data;
 } *subevtdata_ptr;

#endif /* __SUBEVTSTRUCTURES_H__ */
