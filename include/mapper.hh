#ifndef IVANP_MAPPER_HH
#define IVANP_MAPPER_HH

template <typename Cont, typename Cmp>
class mapper {
  mutable Cmp f;
  const Cont* c;
public:
  mapper() = delete;
  template <typename P>
  mapper(const Cont* c, P&& f): f(std::forward<P>(f)), c(c) { }
  mapper(const mapper& m): f(m.f), c(m.c) { }
  mapper(mapper&& m): f(std::move(m.f)), c(m.c) { }
  mapper& operator=(const mapper& m) { f = m.f; c = m.c; return *this; }
  mapper& operator=(mapper&& m) { f = std::move(m.f); c = m.c; return *this; }

  template <typename T>
  auto index(const T& key) const {
    const auto n = c->size();
    std::decay_t<decltype(n)> i = 0;
    for (; i<n; ++i) if (f((*c)[i],key)) break;
    return i;
  }

  template <typename T>
  auto iter(const T& key) const {
    const auto end = c->end();
    auto it = c->begin();
    for (; it!=end; ++it) if (f(*it,key)) break;
    return it;
  }

  template <typename T>
  auto& operator[](const T& key) const {
    for (auto& x : c) if (f(x,key)) return x;
    throw std::out_of_range(cat("nothing mapped to ",key));
  }
};

template <typename Cont, typename Cmp>
inline auto map(const Cont& c, Cmp&& f) {
  return mapper<Cont,Cmp>(&c,std::forward<Cmp>(f));
}

#endif
