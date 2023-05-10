#ifndef MESSAGE_H
#define MESSAGE_H

#include <cfemdef.h>

namespace Message {
  void init      (int* argc, char*** argv, int_t& nproc, int_t& iproc);
  void stop      ();
  void sync      ();
  
  void send      (real_t* data, const int_t N, const int_t tgt);
  void send      (int_t*  data, const int_t N, const int_t tgt);
  void recv      (real_t* data, const int_t N, const int_t src);
  void recv      (int_t*  data, const int_t N, const int_t src);

  void grid      (const int_t& npart2d, int_t& ipart2d,
		  int_t& npartz,        int_t& ipartz);
  
  void exchange  (real_t* data, const int_t nZ,const int_t nP,const int_t sign);
  void exchange  (int_t*  data, const int_t nZ,const int_t nP,const int_t sign);
}
#endif



