#ifndef IVANP_BY_FUNCTION_OR_KEY
#define IVANP_BY_FUNCTION_OR_KEY

#include "type_traits.hh"

namespace ivanp {

template <typename Target, typename Pred, auto_enable_if_t<
  is_indexable<Target,Pred>::value
>...> inline decltype(auto) by_key(Target&& target, Pred&& p)
{ return target[p]; }

template <typename Target, typename Pred, auto_enable_if_t<
  !is_indexable<Target,Pred>::value &&
  has_mf_at<Target,Pred>::value
>...> inline decltype(auto) by_key(Target&& target, Pred&& p)
{ return target.at(p); }

template <typename Target, typename Pred, auto_enable_if_t<
  !is_indexable<Target,Pred>::value &&
  !has_mf_at<Target,Pred>::value &&
  is_callable<Target,Pred>::value
>...> inline decltype(auto) by_key(Target&& target, Pred&& p)
{ return target(p); }

template <typename Target, typename Pred, auto_enable_if_t<
  is_callable<Pred,Target>::value
>...> inline decltype(auto) by_function_or_key(Target&& target, Pred&& p)
{ return p(target); }

template <typename Target, typename Pred, auto_enable_if_t<
  !is_callable<Pred,Target>::value
>...> inline decltype(auto) by_function_or_key(Target&& target, Pred&& p)
{ return by_key(target,p); }

} // end namespace ivanp

#endif
