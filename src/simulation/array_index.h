#ifndef ARRAY_INDEX_H
#define ARRAY_INDEX_H

#include <stdint.h>

static inline 
int32_t _l_array_index(const int32_t xx, 
		       const int32_t yy, 
		       const int32_t zz, 
		       const int32_t ngy,
		       const int32_t ngz) {
  return (xx*ngy + yy)*ngz + zz;
}

static inline 
int32_t _l_inv_array_index_z(const int32_t indx,
			     const int32_t ngy,
			     const int32_t ngz) {
  return indx%ngz;
}

static inline 
int32_t _l_inv_array_index_y(const int32_t indx,
			     const int32_t ngy,
			     const int32_t ngz) {
  return (indx/ngz)%ngy;
}

static inline 
int32_t _l_inv_array_index_x(const int32_t indx,
			     const int32_t ngy,
			     const int32_t ngz) {
  return (indx/ngz)/ngy;
}

static inline int64_t _g_array_index(const int64_t xx, 
				     const int64_t yy, 
				     const int64_t zz, 
				     const int64_t ngy,
				     const int32_t ngz) {
  return (xx*ngy + yy)*ngz + zz;
}

#endif
