#ifndef __FUN4ALLROLLOVERFILEOUTSTREAM_H__
#define __FUN4ALLROLLOVERFILEOUTSTREAM_H__

#include <Fun4AllFileOutStream.h>

class Fun4AllRolloverFileOutStream: public Fun4AllFileOutStream
{
 public:
  Fun4AllRolloverFileOutStream(const std::string &frule = "OUTDATA-%010d-%04d.PRDFF",
			       const unsigned int sizeInMB =0,
			       const int offset =0,
			       const int increment=1,
			       const std::string &name = "Fun4AllRolloverFileOutStream");
  virtual ~Fun4AllRolloverFileOutStream() {}
  int WriteEventOut(Event *evt);
  void identify(std::ostream &os = std::cout) const;

 protected:
  unsigned long long max_file_size;
  int current_sequence;
  int i_offset;
  int i_increment;
  

};

#endif /* __FUN4ALLFILEOUTSTREAM_H__ */
