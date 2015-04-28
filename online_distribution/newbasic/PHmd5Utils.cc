
#include "stdio.h"

#include "PHmd5Utils.h"
#include "md5.h"

int PHmd5File(const char * filename,  unsigned char *digest, int *filesize)
{
  int status;
  FILE *fp;
  fp = fopen (filename, "r");

  if (!fp) return 1; 

  status = PHmd5Stream(fp, digest, filesize);
  fclose(fp);
  return status;
}


int PHmd5Stream(FILE *stream,  unsigned char *digest, int *filesize)
{
  /* Important: BLOCKSIZE must be a multiple of 64.  */
#define BLOCKSIZE 4096
  md5_state_t state;
  char buffer[BLOCKSIZE + 72];
  size_t sum;

  /* Initialize the computation context.  */
  md5_init(&state);

  /* Iterate over full file contents.  */
  *filesize = 0;
  while (1)
    {
      /* We read the file in blocks of BLOCKSIZE bytes.  One call of the
	 computation function processes the whole buffer so that with the
	 next round of the loop another block can be read.  */
      size_t n;
      sum = 0;

      /* Read block.  Take care for partial reads.  */
      do
	{
	  n = fread (buffer + sum, 1, BLOCKSIZE - sum, stream);

	  sum += n;
	  *filesize += n;
	}
      while (sum < BLOCKSIZE && n != 0);
      if (n == 0 && ferror (stream))
        return 1;
      
      /* If end of file is reached, end the loop.  */
      if (n == 0)
	break;
      
      /* Process buffer with BLOCKSIZE bytes.  Note that
	 BLOCKSIZE % 64 == 0
       */
      md5_append(&state, (const md5_byte_t *)buffer,BLOCKSIZE );
	/*      md5_process_block (buffer, BLOCKSIZE, &ctx); */
    }

  /* Add the last bytes if necessary.  */
  if (sum > 0) md5_append(&state, (const md5_byte_t *)buffer,sum );
  /*   md5_process_bytes (buffer, sum, &ctx); */

  /* Construct result in desired memory.  */
  md5_finish(&state,  (md5_byte_t *) digest);
  /*  md5_finish_ctx (&ctx, resblock); */
  return 0;
}

