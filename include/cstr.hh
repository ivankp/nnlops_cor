#ifndef IVANP_CSTR_HH
#define IVANP_CSTR_HH

#include <cstring>

template <size_t N>
bool starts_with(const char* str, const char(&prefix)[N]) {
  for (unsigned i=0; i<N-1; ++i)
    if (str[i]=='\0' || str[i]!=prefix[i]) return false;
  return true;
}

#endif

