#ifndef MAT_HH
#define MAT_HH

template <typename T>
class lt_mat { // lower triangular matrix
  unsigned _n;
  static constexpr unsigned N(unsigned n) { return n*(n-1) >> 1; }
public:
  std::vector<T> mat;
  lt_mat(unsigned n): _n(n), mat(N(_n)) { }
  lt_mat(const lt_mat& o): _n(o._n), mat(o.mat) { }
  lt_mat(lt_mat&& o): _n(o._n), mat(std::move(o.mat)) { }
  T& operator()(unsigned i, unsigned j) {
    if (j>i) std::swap(i,j);
    return mat[N(i)+j];
  }
  const T& operator()(unsigned i, unsigned j) const {
    if (j>i) std::swap(i,j);
    return mat[N(i)+j];
  }
  inline T& operator[](unsigned k) { return mat[k]; }
  inline const T& operator[](unsigned k) const { return mat[k]; }
  constexpr unsigned size() const { return _n; }
};

template <typename T>
class sym_mat { // symmetric matrix
  unsigned _n;
  static constexpr unsigned N(unsigned n) { return n*(n+1) >> 1; }
public:
  std::vector<T> mat;
  template <typename Indexable, typename Pred>
  sym_mat(const Indexable& xx, Pred f): _n(xx.size()), mat(N(_n)) {
    unsigned k = 0;
    for (unsigned i=0; i<_n; ++i) { // (i >= j)
      for (unsigned j=0; j<=i; ++j) {
        mat[k] = f(xx[i])*f(xx[j]);
        ++k;
      }
    }
  }
  sym_mat(const sym_mat& o): _n(o._n), mat(o.mat) { }
  sym_mat(sym_mat&& o): _n(o._n), mat(std::move(o.mat)) { }
  sym_mat& operator+=(const sym_mat& rhs) {
    if (_n!=rhs._n) throw std::length_error(
      "adding sym_mat of unequal size");
    for (unsigned i=mat.size(); i; ) { --i;
      mat[i] += rhs.mat[i];
    }
    return *this;
  }
  T& operator()(unsigned i, unsigned j) {
    if (j>i) std::swap(i,j);
    return mat[N(i)+j];
  }
  const T& operator()(unsigned i, unsigned j) const {
    if (j>i) std::swap(i,j);
    return mat[N(i)+j];
  }
  inline T& operator[](unsigned k) { return mat[k]; }
  inline const T& operator[](unsigned k) const { return mat[k]; }
  constexpr unsigned size() const { return _n; }
};

template <typename T>
struct corr_stdev {
  lt_mat<T> corr;
  std::vector<T> stdev;
};

template <typename T>
corr_stdev<T> corr(const sym_mat<T>& cov) {
  const unsigned n = cov.size();
  lt_mat<T> corr(n);
  std::vector<T> stdev(n);

  for (unsigned i=0; i<n; ++i)
    stdev[i] = std::sqrt(cov(i,i));

  unsigned k = 0;
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j) {
      // std::cout << cov[k] << "/(" << stdev[i]<<'*'<<stdev[j]<<')'<< std::endl;
      corr(i,j) = cov[k]/(stdev[i]*stdev[j]);
      ++k;
    }
    ++k;
  }

  return { std::move(corr), std::move(stdev) };
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const sym_mat<T>& m) {
  for (unsigned i=0, n=m.size(); i<n; ++i) {
    for (unsigned j=0; j<n; ++j) {
      s << std::setw(10) << m(i,j) << ' ';
    }
    s << '\n';
  }
  return s;
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const lt_mat<T>& m) {
  for (unsigned i=0, n=m.size(); i<n; ++i) {
    for (unsigned j=0; j<n; ++j) {
      s << std::setw(10);
      if (j<i) s << m(i,j);
      else s << 0;
    }
    s << '\n';
  }
  return s;
}

#endif
