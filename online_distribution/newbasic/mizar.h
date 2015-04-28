#ifndef __MIZAR_H__
#define __MIZAR_H__
 
typedef struct sdm{
 int sdmlength;
 /* start of mizar data structures */
 unsigned int head0;
 /* conversion 1 information */
 unsigned int head1_tag;
 unsigned int head1;
 
 struct {
     struct {
       unsigned int high_gain[4];
       unsigned int low_gain[4];
     } chip[3];
 } conv1_sm[6];
 /* conversion 2 information */
 unsigned int head2_tag;
 unsigned int head2;
 struct {
    struct {
      unsigned int tac[4];
      unsigned int high_gain[4];
      unsigned int low_gain[4];
    } chip[3];
 } conv2_sm[6];
} *sdm_ptr;
 
typedef struct miz_subdef{
 int sub_length;
 int sub_id;
 int sub_type;
 int sub_decoding;
 /* end of header info*/
 int status;
 int retries;
 int arr[99999];
} *miz_sub;
 
typedef struct miz_irdgdef{
  struct {
   int conv1_high[144];
   int conv1_low[144];
   int conv2_high[144];
   int conv2_low[144];
   int tac[144];
   int trigger[6];
   int write_cell;
   int c1_cell;
   int c2_cell;
   int opcode1;
   int opcode2;
   int dspmap;
   int board_adr;
   int ser_ret;
   int words;
   int byte_err;
   int dummy[14]; /* to round up to 750 words */
  } out[11];
} *miz_irdg;
 
typedef struct miz_indgdef{
 struct {
   int high[144];
   int low[144];
   int tac[144];
   int trigger[6];
   int write_cell;
   int c1_cell;
   int c2_cell;
   int board_adr;
   int ser_ret;
   int socket;
   int port;
   int dsp;
   int words;
   int dummy[3]; /* to round up to 450 words */
 } out[11];
} *miz_indg;
 
typedef struct sdm_c_blockdef{
  int sdmlength;
  unsigned int conv1_info;
  unsigned int conv2_info;
  unsigned int dspmap;
  unsigned int trigger[2];
  int array[999999];
} *sdm_c_block;
 

#endif /* __MIZAR_H__ */
