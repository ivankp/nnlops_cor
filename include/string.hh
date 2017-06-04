#ifndef IVANP_STRING_HH
#define IVANP_STRING_HH

template <typename Str, unsigned N>                                             
bool starts_with(const Str& str, const char(&prefix)[N]) {                      
  for (unsigned i=0; i<N-1; ++i)                                                
    if (str[i]=='\0' || str[i]!=prefix[i]) return false;                        
  return true;                                                                  
}

#endif

