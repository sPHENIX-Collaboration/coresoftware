#ifndef __SIMPLERANDOM_H__
#define __SIMPLERANDOM_H__

/*
 *  After becoming frustrated with the lack of a standalone, portable,
 *  decent random number generator, I decided to make one based on a
 *  cryptographic one-way hash function.  I chose MD5 since it is fast
 *  and free source was readily available.  More cryptographically
 *  secure hash functions are available (e.g. SHA-1), but for the
 *  purposes of a rand/random/erand48 replacement, MD5 should be more
 *  than sufficient.
 *
 *  MD5 takes an arbitrary amount of input and yields a 16 byte hash.
 *  This RNG continually MD5's a 16 byte digest, and uses the bottom N
 *  bits as the random number yielded, where N is just large enough to
 *  include the largest random number desired.
 *
 *      To yield a random number between 0 and r:
 *
 *              create mask which has enough bits to include all of r
 *                      (for example, if r is 100, mask would be 0x7F)
 *
 *              do {
 *                      digest = MD5(digest)
 *                      number = digest & mask
 *              } while (number > r)
 *
 *  The digest should be loaded and saved to a disk file between
 *  invocations of a program using the RNG.
 *
 *  Random functions appear after the included MD5 code.
 *
 *  Send comments to:  skrenta@pbm.com (Rich Skrenta)
 */


/*****************************************************************/

/*
 * This code implements the MD5 message-digest algorithm.
 * The algorithm is due to Ron Rivest.  This code was
 * written by Colin Plumb in 1993, no copyright is claimed.
 * This code is in the public domain; do with it what you wish.
 *
 * Equivalent code is available from RSA Data Security, Inc.
 * This code has been tested against that, and is equivalent,
 * except that you don't need to include two pages of legalese
 * with every copy.
 *
 * To compute the message digest of a chunk of bytes, declare an
 * MD5Context structure, pass it to MD5Init, call MD5Update as
 * needed on buffers full of bytes, and then call MD5Final, which
 * will fill a supplied 16-byte array with the digest.
 */


typedef unsigned int word32;
typedef unsigned char byte;

struct xMD5Context {
        word32 buf[4];
        word32 bytes[2];
        word32 in[16];
};

class simpleRandom {
public:
  simpleRandom();
  simpleRandom(const int iseed);
  
  float gauss(const float mean, const float sigma);
  float rnd(int low, int high);

private:
  
  void byteSwap(word32 *buf, unsigned words);
  void xMD5Init(struct xMD5Context *ctx);
  void xMD5Update(struct xMD5Context *ctx, byte const *buf, int len);
  void xMD5Final(byte digest[16], struct xMD5Context *ctx);
  void xMD5Transform(word32 buf[4], word32 const in[16]);
  void MD5(byte *dest, const byte *orig, int len);
  void load_seed();

  unsigned int digest[4];
    
};
#endif

